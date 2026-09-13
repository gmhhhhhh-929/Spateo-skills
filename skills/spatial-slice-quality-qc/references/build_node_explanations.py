"""Maintain Chinese input/process/output explanations for every workflow node."""
from pathlib import Path
from collections import Counter
import json, xml.etree.ElementTree as ET, html, re
R=Path(__file__).resolve().parent
SVG=R/('spateo_referee_complete_editable.svg' if (R/'spateo_referee_complete_editable.svg').exists() else 'workflow.svg')
OBJECTS=SVG.with_suffix('.objects.json')
nodes={}
def add(id,title,inp,op,out,next,example='',caution=''):
 nodes[id]=dict(id=id,title=title,input=inp,operation=op,output=out,next=next,example=example,caution=caution)

add('a_input','输入数据采用什么组织方式',
 '原始H5AD文件或一个包含多个H5AD的目录；H5AD通常保存每点表达、点的元数据以及空间坐标。',
 '先确定文件里有哪些生物样本、每个点属于哪张切片，再选择相应的切片识别方式。不同样本各自建序列，不把多个发育时期接成一条序列。',
 '待解析的文件列表和数据集边界；此时还没有质量分数。',
 '按数据组织方式进入a_Slice_labels、a_Discrete_z或a_One_file__slice。',
 '7dpa1和10dpa1分别计算，不能让7dpa1最后一张的邻居变成10dpa1第一张。')
add('a_Slice_labels','使用现成的切片标签',
 '一个H5AD及指定的obs切片字段，例如每个点的slice_id。obs是一行对应一个点的元数据表。',
 '按字段值把点分成切片，并按指定顺序或解析出的顺序排列。',
 '每个点的切片归属、唯一切片ID列表和切片顺序。',
 '进入a_contract，核对排序与坐标。','同一文件中slice_id为SL01的所有点组成SL01切片。','标签名称本身不一定表示真实物理顺序。')
add('a_Discrete_z','由离散z值识别切片',
 '带离散第三坐标z的空间数组；每个切片内的z应能形成明确的层。',
 '按z层分组并形成z标签，同时保留实际z和原始x/y。',
 '点到z切片的映射、z切片ID与顺序；x/y继续用于二维QC。',
 '进入a_contract。','z=80这一层可显示为z80；z不进入二维密度面积计算。','连续散布的三维细胞不能未经核对就假定为离散切片。')
add('a_One_file__slice','每张切片一个文件',
 '一组切片H5AD路径，以及可选的path/order清单。',
 '默认按自然文件名顺序排列；若文件名顺序不符合实际切片顺序，使用清单指定。',
 '有序文件/切片列表，每片保留来源文件。',
 '进入a_contract。','自然排序将SL2放在SL10前；普通字典序可能相反。')
add('a_contract','确认切片顺序、原始坐标和来源',
 '上游解析的切片归属、排序、坐标候选键、表达层与元数据。',
 '确认实际选用的原始x/y、顺序和表达层；记录单位信息、点身份和源文件校验值。未知单位应标为未知，不能凭图形大小猜单位。',
 '可追溯的输入契约：每个点是谁、在哪张切片、使用哪套坐标和表达值。',
 '必需字段解析失败→a_stop；通过→a_signals。','spatial与spatial_artificial_unaligned是不同坐标键，不能任意替换。','此节点不执行配准，也不改变原始数据。')
add('a_stop','停止当前受阻的数据集',
 '必需字段、坐标或切片身份无法可靠解析的错误。',
 '记录失败并停止依赖这些信息的当前运行。其他独立数据集是否继续由批处理入口处理。',
 '错误/失败状态；没有可供解释的有效QC分数。',
 '先修正输入契约，再从a_input开始。','无法确定哪列代表切片时，不把每个细胞当作一张切片。')
add('a_signals','检查有哪些可用信号',
 '已确认契约的坐标、表达矩阵和可选细胞类型标签。',
 '分别检查几何、矩阵、注释三类信号。这里是可并行的可用性分流，不是三选一。',
 '可用信号清单及各自处理路径。',
 '坐标→a_geometry；表达矩阵→a_matrix_type；注释→a_labels。','只有可靠坐标也能计算几何证据，但不能凭空生成表达捕获证据。')
add('a_geometry','从坐标计算组织几何',
 '每片原始x/y；若可用，局部低捕获指标还需要该点的捕获计数。',
 '用有限坐标计算凸包面积、点密度、最近邻间距及尾部、连通成分、碎片、空洞等。退化点集按函数规则返回受限值。',
 '每片一组几何指标，而非图像是否“漂亮”的判断。',
 '汇入a_output，随后用于B/C及D中的密度、损伤和连续性域。','最大连通成分占比0.9表示90%的有效几何点属于同一最大连通组织成分。','QC几何中的min_points_geometry=20未实际生效；展示配准中的最小20点是另一条真实规则。')
add('a_labels','从细胞类型标签提取组成信号',
 '真实存在且被选用的细胞类型/注释字段，逐点对应。',
 '统计每片各类型的数量和比例，后续比较相邻切片的组成分布；没有标签则不产生该证据。',
 '每片的类型比例向量；不是基因表达矩阵。',
 '汇入a_output，并参与d_domain_Continuity；标签还可在显式配置时用于展示配准候选排序。','一片神经细胞20%、肌肉细胞80%是一组组成比例。','类型比例改变既可能来自技术缺陷，也可能来自真实解剖变化。')
add('a_matrix_type','判断所选矩阵的语义',
 '选定的表达层/AnnData.X、计数形态检查、元数据及声明的表示类型。',
 '区分计数矩阵、非计数矩阵，以及契约明确声明的注释one-hot X。one-hot即每点以0/1表示所属类别，不是测得的基因计数。',
 '决定允许从矩阵中提取哪些表达/捕获量。',
 '计数→a_semantic_20；非计数→a_semantic_425；声明的one-hot X→a_semantic_830。','每行只有一个1的类型编码，行和总是1，不能据此说每个细胞都捕获了1个转录本。','数值看起来像整数只是线索，不等于生物学计数语义已确认。')
add('a_semantic_20','处理测得的计数矩阵',
 '每行一个点、每列一个基因的非负计数层。',
 '逐点求行和作为捕获计数，求非零基因数；按切片汇总。所选计数层的实际值优先于过时obs汇总。',
 '每点计数/基因数、切片捕获指标及可用于连续性计算的表达汇总。',
 '汇入a_output。','一行[2,0,5]给出总计数7、检出基因2。','捕获计数不等于原始测序reads，也不直接反映测序饱和度。')
add('a_semantic_425','处理归一化等非计数矩阵',
 '非计数表达矩阵，以及可能存在的obs原始计数/基因数汇总。',
 '当前实现优先使用obs计数/基因数；若无obs则使用矩阵行和/非零数作为代理量，并保留语义限制。表达profile仍可进入连续性处理。',
 '注明来源的汇总或代理量；不能无条件解释为真实捕获深度。',
 '汇入a_output。','log-normalized矩阵的行和不是UMI总量；若obs保存了可靠total_counts，应识别并使用该来源。','图画的是当前实际行为；没有表示非计数输入已得到与原始计数等价的验证。')
add('a_semantic_830','处理声明的注释one-hot X',
 '所选矩阵为X，且输入契约明确声明expression_representation包含one-hot。',
 '将捕获、相关表达统计和线粒体代理设为不可用，并跳过表达pseudobulk；保留原始几何和实际注释信息。',
 '几何/注释可用，表达捕获不可用的指标记录。',
 '汇入a_output。','每个点的one-hot类别向量不会被用来推断该点是否低表达。','需要可靠声明；未声明的one-hot矩阵不能假定实现必然识别。')
add('a_output','生成逐片指标和可用性记录',
 '几何、表达和类型分支的结果，以及来源契约。',
 '以切片为单位合并，保留每项指标的原始值、可用性和输入警告。',
 '一张每行对应一片的指标表，以及表达profile/类型向量、点抽样和来源记录等配套对象。',
 '送入B进行邻域比较；每项测量在C转为异常证据。','median_total_counts是一片中所有点总计数的中位数，不是已经标准化的异常分数。')

add('b_window_mode','决定使用固定窗口还是自动选窗',
 '有序切片指标，以及用户配置的window。',
 '固定数字直接指定主窗口宽度；auto比较候选宽度。宽度包含当前片本身。',
 '将采用的选窗方式。',
 '固定→b_fixed；auto→b_auto。','W=3表示前一片、当前片、后一片，共3片。','自动选窗不表示.129/.540/.700这些最终阈值也会按数据集自动改变。')
add('b_fixed','使用固定主窗口',
 '指定的奇数窗口宽度W，默认3。',
 '保留该宽度，不进行候选窗口效用比较；靠近边界时只用实际存在的邻居。',
 '主窗口宽度W及记录。',
 '进入b_sliding_window_example，再按b_context处理每片每指标。','第一张切片的W=3窗口只有自己和后一片，不虚构前一片。')
add('b_auto','自动比较候选窗口',
 '整组原始指标和候选3/5/7，候选必须不超过系列长度。',
 '对每个候选宽度，最多选择8个确定性位置施加有限的表达指标扰动，计算检出、得分增量和基线非keep比例，并对较大窗口加惩罚。按记录的效用选择。',
 '一个主窗口宽度及各候选效用记录。没有可用候选时回退W=3并注明系列过短。',
 '进入b_sliding_window_example。','效用=扰动检出比例+0.35×平均分数增量−0.35×基线非keep比例−0.012×(W−3)。','这是选窗启发式，不是独立准确率评估；这里最多8个位置不是8组数据或8轮正式验证。')
add('b_sliding_window_example','滑动窗口示意中的切片、虚线框与箭头',
 '按切片索引排序的系列，以及主窗口W。',
 '紫色i是当前片，虚线框是当前窗口；计算参照时排除i，只取框内有效邻居。随后i前进一位，重新为下一张切片计算。',
 '每张切片自己的邻居集合；没有一个固定的全局“正常切片”充当所有参照。',
 '进入b_context。','图中W=5时，i的候选邻居是i−2、i−1、i+1、i+2；最外侧省略号不是具体数据点。','这里是切片顺序窗口，不是空间坐标里的ROI框，也不是跨物种窗口。')
add('b_context','逐指标检查有效邻居位于哪一侧',
 '当前片、候选邻居，以及某一指标的有限数值。',
 '分别查看当前片前后是否有有效值。同一切片不同指标可能因缺失情况不同而走不同分支。',
 '双侧、仅单侧或无有效邻居的状态。',
 '双侧→b_context_25；单侧→b_context_425；无邻居→b_context_825。','两侧切片都存在，但前一片的计数指标缺失时，该计数指标可能只能使用后一侧。')
add('b_context_25','用双侧邻居估计局部期望',
 '当前片前后都有该指标的有效值。',
 '排除当前片，对邻居数值按位置拟合中位两两斜率趋势，再在当前索引处求期望；同时记录左右中位数差异。',
 '局部期望e、双侧标记及左右不一致程度。',
 '在b_support安排主/支持窗口计算，再送C标准化。','若前后正常计数分别约900和1100，当前片的期望可能约1000；不把当前片自身低计数混入参照。')
add('b_context_425','用单侧邻居估计期望',
 '端点或缺失后，只剩当前片一侧的有效邻居。',
 '使用这些邻居的中位数；窗口在边界截断，不循环到序列另一端，也不补造切片。',
 '单侧期望e和单侧标记；后续指标证据与总分受到边界减权。',
 '继续b_support/C；最终E不能仅靠单侧证据自动排除。','第一张切片只有后方邻居；不是因为缺少前方切片就认定它质量差。')
add('b_context_825','没有可用参照',
 '窗口内该指标没有任何有效邻居。',
 '期望保留为不可用，该指标内部异常值为0；不补值，也不制造可比较性。',
 '无参照的记录和“未获得异常证据”的内部值。',
 '继续合并其他可用指标；该项不能作为排除依据。','所有邻居都没有表达测量时，不能用“表达异常=0”推断当前片表达正常。')
add('b_support','合并主窗口与稳健支持窗口的证据',
 '主窗口结果，以及同组数据可支持的较宽趋势窗口。',
 '分别重算指标异常，并保留两者较大值。常用配置下支持宽度最多9；实际公式会保证不小于主窗口。主窗口的双侧信息继续用于切片上下文。',
 '最终单指标异常、两套期望/证据记录，以及切片上下文与score_confidence。',
 'C描述这里使用的单指标算法；合并后的证据进入D。','连续几张低质量片可能互相掩盖；较宽支持趋势提供另一组参照。','B→C是说明分工：实现会在主窗口和支持窗口分别调用C，再取max，不是只计算一次C。切片双侧字段为各评分指标双侧状态取any；confidence为1或.72，不是正确概率。')

add('c_input','准备当前观测与局部期望',
 '原始指标值v，以及B估计的同一指标期望e。',
 '确保比较的是同一测量、同一单位；检查有限值。',
 '待计算的差值、方向和可比较状态。',
 '进入c_direction。','v=400、e=1000可表示当前与期望的每点计数中位数；它们都还不是异常分数。')
add('c_direction','选择异常的方向',
 '当前指标名称及预定义的下降型/上升型规则。',
 '按照指标含义选择“减少不利”还是“增加不利”；方向不是看这次数据升降后临时决定。',
 '相应的方向性效应计算方式。',
 '下降型→c_low；上升型→c_high。','计数下降可提示捕获损失；空洞比例上升可提示组织内部缺失。')
add('c_low','计算下降型相对效应',
 '下降型指标的观测v和期望e。',
 'm=max((e−v)/|e|,0)，分母用极小正数保护。只给不利下降记正效应。',
 '相对下降幅度m，以及供稳健残差使用的e−v。',
 '同时进入c_effect与c_z。','e=1000、v=400，则m=.60；v=1200时这一路m=0。','这不是把全部指标都除以最大值进行min-max归一化。')
add('c_high','计算上升型绝对效应',
 '上升型指标的观测v和期望e。',
 'm=max(v−e,0)。比例型指标接近0时，使用绝对增量避免相对倍数爆炸。',
 '绝对增加量m，以及供稳健残差使用的v−e。',
 '同时进入c_effect与c_z。','空洞比例由.02升为.08，m=.06（6个百分点），不是异常分数等于3或4。')
add('c_effect','把效应大小映射到0–1',
 '效应m，以及该指标预设的mild、severe两个效应界限。',
 'R=clip((m−mild)/(severe−mild),0,1)。clip表示小于0设0、大于1设1。',
 '效应证据R；并非最终总分。',
 '与c_z的结果在c_merge合并。','计数的mild=.25、severe=.62；相对下降.60对应R约.946。','这里mild/severe是单指标效应界限，与第一阶段.129/.700及域阈值.65/.60不同。')
add('c_z','计算稳健残差异常',
 '整个序列中该指标的方向性残差，以及当前片残差。',
 '用残差中位数居中，以1.4826×MAD估计尺度，退化时用标准差，再退化为1。正向标准化残差z映射为Z=clip((z−1.5)/3,0,1)。',
 '当前变化相对整组残差分布有多突出，形成证据Z。',
 '上升型先经c_small_high；下降型直接进入c_merge。','两个指标有相同绝对变化，但平时波动尺度不同，稳健残差证据可以不同。','MAD是“绝对偏差的中位数”；这不是检验p值，也不是以当前片自己的点分布估计缺陷概率。')
add('c_small_high','限制近零指标被夸大的异常',
 '上升型的Z、绝对效应m及轻效应界限mild。',
 '把Z乘clip(m/mild,0,1)。当绝对变化很小时，即使残差标准化值大，也不给它满额证据。',
 '经过最小效应约束的上升型Z。',
 '进入c_merge；下降型旁路跳过本节点。','m=.004、mild=.04，则最多保留原Z的10%。','只作用于上升型稳健残差分支，不同时缩小R。')
add('c_merge','合并证据并降低不可靠邻域的权重',
 'R、处理后的Z、双侧/单侧状态，以及左右邻居的不一致程度。',
 '先取max(R,Z)，再乘上下文因子（双侧1、单侧.72）和邻居一致性惩罚（.35–1）。计数/基因另外保留直接相对下降证据；主/支持窗口结果再取max。',
 '每项指标的异常A及其期望、支持和缺失信息，供四域聚合使用。',
 '跨区进入d_domain_Density、d_domain_Expression、d_domain_Damage和d_domain_Continuity。','若R=.4、Z=.8、双侧且一致性惩罚=.5，则这一路A=.4；计数直接效应或较宽窗口可能再提高最终单指标证据。','表达profile与类型组成有专用散度映射，不能直接把本节点这条通用式套到所有连续性信号。')

add('d_domain_Density','密度与组织量证据域D',
 '来自C的密度、点数、空洞、近邻尾部、碎片及最大成分占比异常分数。',
 '按图中权重聚合，并与最强两个指标的组合取较大值。这样既看总体，也避免少数强信号被均值稀释。',
 '0–1的密度/组织量域分数D。',
 '通过d_available所示共同规则，进入d_rawscore。','输入的“密度异常=.8”已是标准化分数，不再是每平方微米的细胞数。','这一域共享了碎片、空洞等信号，因此不能把它和损伤/连续性视作完全独立的证据。')
add('d_domain_Expression','表达捕获证据域X',
 '计数、基因、零值比例、复杂度、检出基因比例、局部/区域低捕获的异常分数。',
 '计算可用指标的加权均值，与最强两个指标组合取较大值。',
 '0–1的表达捕获域分数X。',
 '进入d_rawscore及解剖保护/检测器规则。','计数和检出基因同时明显降低，可以共同支持这个域。','图中的X是一个标量域分数，不是原始表达矩阵AnnData.X。')
add('d_domain_Damage','损伤代理证据域M',
 '线粒体比例异常和组织碎片异常。',
 '对可用异常分别赋权.65、.35并聚合。',
 '0–1的损伤代理域分数M。',
 '进入d_rawscore和解剖保护规则。','线粒体比例异常可能支持应激/损伤，但单靠这个指标不能确认具体实验原因。','这是损伤代理，不是组织损伤的直接诊断。')
add('d_domain_Continuity','跨切片连续性证据域C',
 '表达profile异常、细胞类型组成异常，以及空洞异常。',
 '表达profile先比较当前片与邻居平均profile的余弦差异，并考虑左右邻居相似度；类型组成用Jensen–Shannon散度。映射成异常后按.58/.22/.20聚合。',
 '0–1的连续性域分数C。',
 '进入d_rawscore和解剖保护规则。','若当前片表达概貌突变，而前后两侧彼此相似，该片可能获得连续性异常证据。','profile是按片汇总的基因表达概貌；此处并未用展示预配准坐标计算逐点解剖对应。')
add('d_available','四域共享的域内聚合规则',
 '某一域中的指标异常分数及固定权重。',
 '可用有限值做加权均值并按可用权重归一；密度/表达域另算.65×最强+.35×次强，再与均值取max。',
 '各域分数的计算规则，不是额外第五个评分域。',
 '得到四个域分数后进入d_rawscore。','最强=.8、次强=.6时，top-two=.73。若该域加权均值=.5，则密度/表达域采用.73。','上游有些不可用信息被内部记为0证据；0不等于观察到正常。“finite可用值”是实现规则，不能据此补出生物学测量。')
add('d_rawscore','从四个域得到未保护总分S₀',
 'D、X、M、C四个域分数。',
 '计算.31D+.39X+.12M+.18C，同时计算最大域支路clip(.48+.78×(最大域−.68),0,1)，两者取max；统计有多少域≥.45。',
 '未施加解剖保护和端点减权的总分S₀、支持域数及各域值。',
 '进入d_protection；总分还不是最终公开动作。','D=.65、X=.55、M=.10、C=.20时，加权值=.464，最大域支路=.4566，故S₀=.464。','S₀不是四个域简单平均，也不是概率；最大域支路解释了为何只降低review分界可能并无效果。')
add('d_protection','是否触发解剖保护',
 '四域分数及点数、密度、空洞的单指标异常。',
 '检查d_guard_rule列出的两种组合，任一成立就设置保护标记。此处想减少把自然收窄或局部结构差异误当技术缺陷的风险。',
 '保护标记true/false。',
 'yes→d_guard_yes；no→直接将S₀送d_endpoint。','组织量减少，但密度、表达和空洞仍相对正常时，可能符合自然收窄保护。','保护是经验规则，不是已经确认该组织变化一定合理。')
add('d_guard_rule','解剖保护条件定义框',
 '各域和单指标异常分数，不是原始点数或原始密度。',
 '条件一：D≥.48、X<.32、M<.38、C<.45。条件二：点数异常≥.45、密度异常<.30、空洞异常<.30、X<.35。任一组完整成立即保护。',
 '提供给d_protection的逻辑判断定义。',
 '虚线指向d_protection；这不是另一条会修改分数的独立处理步骤。','“点数异常≥.45”不等于“点数≥.45”或“减少45%”，而是标准化后的点数异常分数。')
add('d_guard_yes','对受保护切片限制自动排除',
 '保护标记为true，以及S₀。',
 'S取S₀和.62的较小值；.62来自内部检测器默认排除阈值.64减.02。保留保护标记，后续排除门控也检查它。',
 '保护后的S和标记。',
 '进入d_endpoint。','S₀=.80时被限制为.62；S₀=.40时仍为.40。','保护不是将切片修好；只是限制当前自动判定。')
add('d_endpoint','施加末端减权并交出最终异常分数',
 '保护处理后的S或未保护的S₀，以及主窗口上下文。',
 '单侧上下文总分再乘.82；双侧不乘。最终截断到0–1。',
 '供最终策略使用的异常总分S，同时保留域分数、保护标记、上下文和confidence。',
 '送入e_detector，再由e_stage1决策。','单侧S=.50变为.41；同时仍有单侧标记，不能因为分数够高就绕过双侧排除要求。','单指标阶段的.72和这里总分阶段的.82是两处不同操作。')

add('e_detector','生成内部检测器候选标签',
 'S、各域、单指标异常，以及双侧/保护状态。',
 '按框内条件生成内部keep/review/exclude：强表达或特定几何组合可触发候选；若缺双侧或有保护，不能从这里直接判内部exclude。',
 'recommendation内部标签及原因，不是最终公开final_call。',
 '与同一总分S一并送入e_stage1。','内部阈值.38/.64与最终策略.129/.700属于不同层，不能拿其中一组替代另一组。','此处review意味着存在值得后续处理的证据，不表示最终报告要求用户再做人工裁决。')
add('e_stage1','第一阶段按锁定总分阈值分流',
 '最终异常总分S、内部检测器标签与保护/上下文。',
 'S≤.129走保留候选；.129<S<.700进review；S≥.700还要经过直接排除门控。',
 '阈值区间threshold_band及第一阶段路径。',
 '低→e_keep_candidate；中→e_middle；高→e_high。','S=.50和S=.60都先进入review，随后才按.540区别细筛门槛。','.540不是第一阶段直接keep阈值；当前策略也不按数据集重估这三条边界。')
add('e_keep_candidate','低总分直接保留候选',
 'S≤.129的切片及其上下文记录。',
 '在当前完整二元实验应用中给出保留动作，并保留内部检测器/保护等审计字段；不做第二阶段review证据门控。',
 '最终keep及保留理由。',
 '沿左侧线进入e_keep，再写f_audit。','S=.10不会为了形式完整再走一次低分review。','这里keep不表示“独立认证健康”；当前实验作用域不生成生物学认证。')
add('e_middle','建立内部review队列',
 '.129<S<.700的切片；此外也接收高分直接排除门控未通过的切片。',
 '收集仍需要多维/跨窗口证据判断的切片。',
 '待第二阶段处理的review记录。',
 '进入e_review_split。','S=.45将实际接受低分细筛，而不是一律keep-only提前返回。')
add('e_high','高总分的直接排除门控',
 'S≥.700的切片。',
 '同时要求内部检测器为exclude、主窗口双侧、没有解剖保护。任一失败即降入review。',
 '直接exclude，或降级review及失败原因。',
 '全部通过→e_exclude；任一失败→e_middle/e_review_split。','即使S很高，只要缺乏双侧支持也不会直接排除。','图右侧“all pass”长线专属这条直接排除路径，不表示跳过了保护检查。')
add('e_review_split','为每一条review选择细筛分段',
 '已经进入review的切片及其主窗口总分S。',
 'S<.540进入低分细筛；S≥.540进入高分细筛。这里只选用规则，没有重新计算总分。',
 '当前切片使用的tier名称和参数。',
 'yes→e_low_tier；no→e_high_tier。','S恰好=.540时进入高分段。','这是主窗口分段；确认窗口不再用.540把当前切片切换到另一套规则。')
add('e_low_tier','低分review的证据要求',
 '.129<S<.540的review切片、四域证据和各确认窗口结果。',
 '选用至少2个支持域、最强域≥.65的规则；每个支持域自身必须≥.45。然后还要经过共用门控。',
 '低分段要求，不是仅凭最强域就给exclude。',
 '进入e_all_gates。','最大域=.64时，尽管可能有两个支持域，仍未达到本段最强域门槛。','低分门槛.65高于高分段.60：对这项证据而言，低分需要更强理由才能排除。')
add('e_high_tier','高分review的证据要求',
 'S≥.540且仍处于review的切片，包括直接排除失败后降级的高分片。',
 '选用至少2个支持域、最强域≥.60的规则；继续检查所有共用条件。',
 '高分段规则。',
 '进入e_all_gates。','S=.60、最强域=.62，只通过了分数与最强域这两项，仍可能因某窗口不支持而保留。','本段没有.700上界，因为还包括降级进入review的高分片。')
add('e_all_gates','执行第二阶段全部AND条件',
 '被选中的tier；主窗口域证据、检测器和保护；可用3/5/7确认窗口的完整详情。',
 '主窗口必须是检测器review/exclude、双侧、无保护、confidence≥.90，并满足所选域数/最强域。每个可用窗口还必须满足检测器、双侧、无保护、域数/最强域及该tier总分下限。要求支持比例100%；缺窗口详情视作不通过。',
 '通过/失败、逐条失败理由、窗口支持比例。',
 '全部通过→e_exclude；任一失败→e_keep。','高分段主窗口通过，但7片窗口最强域只有.58时，不能达到.60门槛，最终keep。','确认窗口最低总分为低段.129、高段.540；只主窗口另查confidence。窗口检查不是把图B的auto选窗重新做一遍，也不是按“多数票”通过。')
add('e_keep','最终保留动作',
 '第一阶段保留候选，或review细筛任一必要条件未通过的切片。',
 '写final_call=keep，记录来自哪条路径和为什么不能自动排除。',
 '公开keep标签和完整内部审计。',
 '进入f_audit与viewer。','“跨窗口不稳定而保留”和“总分低而保留”最终同是keep，但审计理由不同。','keep代表当前证据不足以支持排除，不保证不存在缺陷；原始点不会被删除。')
add('e_exclude','最终排除建议',
 '第一阶段直接门控全部通过，或第二阶段全部证据条件通过的切片。',
 '写final_call=exclude，保存直接/细筛路径及证据。',
 '公开exclude标签，供下游挑选输入时参考。',
 '进入f_audit与viewer。','同为exclude，可分别来自高总分直接路径或低分多域稳定细筛路径。','不会自动删除源文件；当前实验策略的exclude也不等于病理或技术原因已经确诊。')

add('f_audit','写出可追溯的逐片决策记录',
 '每片原始指标、S、内部标签、两阶段路径、窗口证据及final_call。',
 '一片一行保存，记录policy和输入校验值。实验作用域保留实际二元动作，但清除独立认证结论。',
 'binary audit/calls CSV、策略应用JSON及后续报告需要的信息。',
 '输出到报告；按需求进入f_prereg_request的独立展示路径。','用户看到keep/exclude；研究者仍可从audit查到review分段和失败条件。','阈值分流、内部候选、最终动作、认证状态是不同字段。')
add('f_prereg_request','是否需要共同坐标系中的比较',
 '显示需求、原始坐标和已有QC结果。',
 '若需要跨片共享ROI，建立或复用展示预配准；否则直接用原始坐标展示。',
 '展示坐标来源的选择。',
 'no→f_raw_only；yes→f_rigid。','只看每片全局统计时可不做预配准；要用同一个框比较邻片时需要明确共同frame。','这条路径不回写D/E中的评分或决策。')
add('f_raw_only','使用原始坐标展示',
 '原始坐标及原QC结果。',
 '直接绘制原始点，不添加共同参考变换。',
 '原始坐标下的证据视图。',
 '进入f_viewer；不据此声称可做共同解剖位置的ROI比较。','独立旋转过的两片即使坐标数值范围相同，也不意味着图中同一位置是同一组织。')
add('f_rigid','轻量化刚性预配准',
 '同一个数据集内有序原始x/y、固定配置、可选显式注释排名条件。',
 '每片确定性抽样最多1000点，24个角度起始，每个最多50迭代，以较近80%对应点拟合平移/旋转。相邻变换正确组合到该数据集中间片参考系。',
 '每片原始→展示矩阵T、参考片、逐对关系、抽样和状态；将T应用到全部原始点产生独立展示坐标。',
 '进入f_fit_status。','若B已映到参考A，C先映到B，则C到A用T(B→A)×H(C→B)，不能只显示C对B的局部坐标。','不缩放、不镜像、不非刚性拉伸；拟合最少20点与A区QC几何配置是不同参数。已有缓存可复用，不必每次重拟合。')
add('f_fit_status','检查拟合失败、方向歧义和参考链状态',
 '拟合矩阵、几何残差、替代角度结果、可选注释支持和上游参考链。',
 '检查是否报错、不同方向解是否接近、几何/注释差异是否过大，以及父级链是否已有警告。当前启发式如替代方向代价比<1.05或归一化几何残差>.08可触发警告。',
 'estimated、warning或failed等状态，以及警告原因和受影响切片。',
 '有失败/歧义→f_warn；获得拟合且无记录警告→f_fit_ok。','两个明显不同旋转方向的拟合代价很接近时，要提示方向歧义。','图上“reliable enough”是状态/启发式检查的概括，不是经过真值验证的解剖准确率。')
add('f_warn','把不可靠预配准显式显示出来',
 '配准失败、歧义或受警告参考链影响的切片。',
 '保留失败和传播状态；显示能提供的坐标及警告，说明局部比较受限，不静默改为成功。',
 '带可见限制的展示状态。',
 '进入f_viewer；用户解读ROI时必须同时看到警告。','若相邻配准边不可靠，经该边组合到参考系的后续片也应保留相关警告。','配准误差不直接变成切片QC的低质量证据。')
add('f_fit_ok','保存已估计的展示坐标',
 '合法刚性矩阵、有效参考链且没有触发当前警告的拟合。',
 '保存矩阵和逆矩阵、展示坐标、身份映射、参数和哈希。',
 '共同frame中的点及原始对应关系。',
 '进入f_viewer。','原始点[10,20]经T变成展示点[12,25]，仍保留原始点身份与原始坐标。','estimated只是算法完成估计；无警告不证明解剖对应完全正确。')
add('f_viewer','生成现有多面板证据报告',
 '二元决策、原始指标、点级表达、KDE/成分证据，以及可选的展示坐标和状态。',
 '组织Overview、Statistics、Slice evidence、All slices、Methods页面，联动3/5片窗口及KDE、表达捕获、成分图。',
 '可检查的HTML报告或多数据集集合；每片显示QC结论、证据与坐标来源。',
 '共享坐标可用时进入f_roi；原始视图仍可独立检查。','KDE图颜色表示局部点密度估计，表达图颜色表示捕获/基因信号；二者不是最终总分热图。','显示3/5片的窗口不等于QC主窗口或3/5/7确认窗口；颜色和视图缩放不修改final_call。')
add('f_roi','准确框选并追溯原始点',
 '指定frame_id中的展示矩形、邻片的原始→展示矩阵和点身份。',
 '将点映到展示frame后做矩形包含判断；把矩形四角经逆变换映回原始坐标，得到多边形。浏览器导出抽样选点，服务器export-roi可校验哈希后在全点缓存重放。',
 '每片选中的原始行号/obs身份、坐标、计数、frame_id和原始ROI多边形。',
 '交付可追溯选区；它不回到D/E重新改变QC。','旋转后的展示矩形在原始坐标里是斜四边形，不能用其轴对齐外接框代替精确选区。','相同frame范围不保证完全相同的解剖区域；抽样点数量也不等于全量ROI点数量。')

# Match explanations to real editable SVG objects, not an independently redrawn figure.
objects=json.loads(OBJECTS.read_text())
objects.insert(next(i for i,o in enumerate(objects) if o['id']=='b_context'),dict(panel='B',id='b_sliding_window_example',type='illustration',x=115,y=540,w=1000,h=200,source='_local_expected; _choose_window'))
assert set(nodes)=={o['id'] for o in objects}, (set(nodes)-{o['id'] for o in objects},{o['id'] for o in objects}-set(nodes))
svg_text=SVG.read_text();root=ET.fromstring(svg_text);ns={'s':'http://www.w3.org/2000/svg'}
groups={g.get('id'):g for g in root.findall('.//s:g',ns)}
counts=Counter();rows=[]
for obj in objects:
 counts[obj['panel']]+=1;code=f"{obj['panel']}{counts[obj['panel']]:02d}"
 row={**obj,**nodes[obj['id']],'code':code}
 g=groups[obj['id']];text=g.find('s:text',ns);row['label']=text.text if text is not None else row['title']
 if not row['source']:row['source']={'A':'_extract_one_file','B':'_choose_window','C':'_metric_anomaly','D':'_score_metrics','E':'apply_high_confidence_policy; _evaluate_review_tier_row','F':'slice_preregistration.py; slice_quality_visualization.py'}[row['panel']]
 rows.append(row)
by_id={r['id']:r for r in rows}
# Replace identifier mentions with human-facing code + short Chinese names.
for row in rows:
 for key in ['next','operation','caution','output']:
  for id in sorted(by_id,key=len,reverse=True):
   row[key]=row[key].replace(id,by_id[id]['code']+'「'+by_id[id]['title']+'」')
(R/'workflow_node_guide.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')

# Add small live-text ID badges to a separate annotated copy, preserving original SVG.
annotated=svg_text
for row in rows:
 id=row['id'];badge=f'''<g id="number_{row['code']}" pointer-events="none"><rect x="{row['x']+7}" y="{row['y']-11}" width="53" height="26" rx="5" fill="#183047"/><text x="{row['x']+33.5}" y="{row['y']+8}" font-family="Arial" font-size="17" text-anchor="middle" font-weight="700" fill="white">{row['code']}</text></g>'''
 # Place the badge after the node geometry so it remains visible and editable.
 start=annotated.index('<g id="'+id+'"');opening_end=annotated.index('>',start)+1
 depth=1
 for match in re.finditer(r'</?g\b[^>]*>',annotated[opening_end:]):
  depth += -1 if match.group().startswith('</') else 1
  if depth==0:
   end=opening_end+match.start();break
 annotated=annotated[:end]+badge+annotated[end:]
annotated=annotated.replace('</svg>','</svg>')
(R/'spateo_referee_annotated.svg').write_text(annotated)

intro='''# Spateo Referee：逐节点输入、操作与输出说明

本说明覆盖完整流程图的60个处理/判断节点，以及1个滑动窗口示意控件，共61项。编号与[带编号SVG](spateo_referee_annotated.svg)及[可点击对照页面](workflow_explained.html)一致；原始未编号SVG保留不动。

## 先明确图中流动的是什么

A输出原始测量和输入来源；B为每片确定比较对象；C把每个测量的偏离程度转为异常证据；D聚合为四域和总分；E产生最终keep/exclude；F把证据、坐标和决策组织为报告。大区箭头表示这些信息的依赖，不表示每一区只执行一次：C会被B的主/支持窗口调用，第二阶段确认窗口也会重复相应评分。

图中的矩形一般表示处理或数据结果，菱形表示选择/判断；A的可用信号是并行分流，非互斥三选一。实线表示信息流或判断后的去向；带yes/no、all pass/any fails的线表示条件成立/失败；D的虚线是规则定义指向判断，不是额外处理。跨区C/D圆标只是续接位置，没有计算功能。浅绿色keep与红色exclude才是公开动作，其余颜色帮助分组，不是独立阈值。

## 常用词与符号

| 名称 | 此处含义 |
|---|---|
| 点 / spot / cell | 数据的一行空间观测，可能是细胞或测量spot，依数据平台而定 |
| H5AD / obs / obsm | 数据文件 / 逐点元数据表 / 坐标等逐点数组 |
| i / W | 当前切片索引 / 含当前片在内的窗口总宽度 |
| v / e | 当前片的原始指标值 / 邻居估计的同一指标局部期望 |
| m / mild / severe | 不利变化效应量 / 单指标轻、重效应界限 |
| R / z / Z / A | 效应证据 / 稳健标准化残差 / 映射到0–1的残差证据 / 合并后的单指标异常 |
| D / X / M / C | 密度、表达、损伤代理、连续性四个域分数；域X不是AnnData.X |
| S₀ / S | 未保护的综合分数 / 保护与端点处理后的最终异常总分 |
| supporting domain | 域分数≥.45，属于支持异常的域；不意味着统计独立 |
| clip(x,0,1) | 将x限制在0到1之间；低于0设0，高于1设1 |
| MAD / SD | 绝对偏差中位数 / 标准差，用于刻画波动尺度 |
| profile / pseudobulk | 将一个切片内点的基因表达汇总形成表达概貌 |
| KDE / NN | 核密度估计 / 最近邻；是空间测量或展示，不是最终质量概率 |
| confidence | 当前邻域支持字段，主要取1或.72，不是正确率概率 |
| tier / AND | 某分数区间采用的一套门控 / 列出的条件必须全部成立 |
| T / H / frame_id | 原始到共同参考的变换 / 相邻片变换 / 展示坐标系唯一标识 |

同样的0.65可能出现在不同公式中：域内top-two的0.65是权重，低分细筛的0.65是最强域门槛。它们不能混用。最终分流.129/.540/.700也不是单指标mild/severe。

## 一条示例贯穿流程（仅用于说明，并非新增实验）

假设一张内部切片有可靠计数和双侧邻居。某指标期望计数1000、当前400，则相对下降m=.60；计数效应映射R约.946，但这仍要结合残差、上下文和其他证据，不能直接当最终总分。若最终四域为D=.65、X=.65、M=.10、C=.20，则S₀=.503；没有保护且双侧时S=.503。第一阶段进review，第二阶段选低分tier，要求2个域≥.45、最强域≥.65及所有其他门控。主窗口和全部可用确认窗口均通过才exclude；任一必要窗口失败则keep。示例没有指定真实切片，也不构成性能结果。

'''
md=[intro];panel_names={'A':'输入与信号选择','B':'窗口与边界','C':'单指标标准化','D':'域聚合与保护','E':'两阶段筛选','F':'审计、展示与ROI'}
last=None
for row in rows:
 if row['panel']!=last:md.append(f"\n## {row['panel']}区：{panel_names[row['panel']]}\n");last=row['panel']
 md.append(f"\n### {row['code']} · {row['title']}\n\n图中文字：**{row['label']}**\n\n")
 for key,label in [('input','输入'),('operation','具体操作'),('output','输出'),('next','后续去向'),('example','例子'),('caution','容易误解的地方')]:
  if row[key]:md.append(f"**{label}：**{row[key]}\n\n")
 md.append(f"实现追溯：`{row['source']}`；SVG节点ID：`{row['id']}`。\n")
(R/'NODE_EXPLANATIONS.md').write_text(''.join(md))

# Self-contained interactive explanation. No external library or browser network request.
svg_inline=annotated[annotated.index('<svg'):]
for row in rows:svg_inline=svg_inline.replace('<g id="'+row['id']+'"', '<g data-guide="'+row['id']+'" tabindex="0" role="button" aria-label="'+html.escape(row['code']+' '+row['title'],quote=True)+'" id="'+row['id']+'"')
encoded=json.dumps(rows,ensure_ascii=False).replace('</','<\\/')
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Spateo Referee · 逐节点说明</title><style>
*{box-sizing:border-box}body{margin:0;font-family:Arial,"PingFang SC",sans-serif;color:#183047;background:#f2f5f8}header{padding:18px 24px;background:white;border-bottom:1px solid #cad5e0}h1{font-size:23px;margin:0 0 8px}header p{margin:7px 0;font-size:14px;color:#53697d}a{color:#245dab}nav{display:flex;gap:20px;font-size:14px}.layout{display:grid;grid-template-columns:minmax(0,1fr) 455px;height:calc(100vh - 130px)}.canvas{overflow:auto;padding:18px;position:relative}.canvas svg{display:block;width:1800px;height:auto;background:white}.tools{position:sticky;top:0;left:0;z-index:2;display:flex;gap:8px;width:max-content;padding:8px;background:#fff;border:1px solid #ccd6df;border-radius:8px}button,select,input{font:inherit;padding:8px;border:1px solid #b9cad8;border-radius:6px;background:white;color:#183047}button{cursor:pointer}.sidebar{background:white;border-left:1px solid #ccd6df;overflow:auto;padding:22px}.sidebar h2{font-size:23px;line-height:1.5;margin:15px 0 4px}.english{font-size:13px;color:#617489;overflow-wrap:anywhere}.sidebar label{font-weight:700;display:block;margin-top:18px}.sidebar p{font-size:15px;line-height:1.8;margin:6px 0;overflow-wrap:anywhere}.sidebar select,.sidebar input{width:100%;margin-bottom:8px}.source{font-size:12px!important;background:#f3f6f9;padding:10px;border-radius:6px}.warning{border-left:3px solid #dca541;padding-left:12px}.route button{margin:5px 5px 0 0;font-size:13px}.pager{display:flex;justify-content:space-between;margin-top:20px}g[data-guide]{cursor:pointer}g[data-guide]:hover>rect,g[data-guide]:hover>path{stroke:#2575a9;stroke-width:4}g[data-guide].selected>rect,g[data-guide].selected>path{stroke:#ce7228;stroke-width:5}g[data-guide]:focus{outline:none}g[data-guide]:focus>rect,g[data-guide]:focus>path{stroke:#2575a9;stroke-width:5}@media(max-width:950px){.layout{grid-template-columns:1fr;height:auto}.canvas{height:55vh}.sidebar{height:45vh}.canvas svg{width:1700px}}@media print{header,.canvas,.tools,.pager,input,select{display:none}.layout{display:block;height:auto}.sidebar{overflow:visible;border:0}}
</style><header><h1>流程图逐节点说明</h1><p>点击任一节点查看“输入 → 操作 → 输出 → 下一步”。共61项。图中编号与文字版一致，原始QC参数和结果未修改。</p><nav><a href="spateo_referee_annotated.svg" download>带编号可编辑SVG</a><a href="NODE_EXPLANATIONS.md">完整文字说明与符号表</a><a href="index.html">返回整图</a></nav></header><div class="layout"><div class="canvas" id="canvas"><div class="tools"><button id="smaller">缩小</button><button id="larger">放大</button><button id="fit">适合宽度</button><span>橙色描边＝当前节点</span></div>__SVG__</div><aside class="sidebar"><input id="search" placeholder="搜索编号、名称、输入或输出" aria-label="搜索节点"><select id="nodeSelect" aria-label="选择节点"></select><div id="detail" aria-live="polite"></div><div class="pager"><button id="previous">上一个</button><button id="locate">定位图中节点</button><button id="next">下一个</button></div></aside></div><script id="guide-data" type="application/json">__DATA__</script><script>
const rows=JSON.parse(document.getElementById('guide-data').textContent), byId=Object.fromEntries(rows.map(r=>[r.id,r]));
const select=document.getElementById('nodeSelect'),detail=document.getElementById('detail');let active=rows[0].id;
function esc(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));}
function populate(term=''){const filtered=rows.filter(r=>(r.code+' '+r.title+' '+r.label+' '+r.input+' '+r.output).toLowerCase().includes(term.toLowerCase()));select.innerHTML=filtered.map(r=>`<option value="${r.id}">${r.code} · ${esc(r.title)}</option>`).join('');if(filtered.some(r=>r.id===active))select.value=active;return filtered.length;}
function show(id){if(!byId[id])return;active=id;const r=byId[id];document.querySelectorAll('g[data-guide]').forEach(g=>g.classList.toggle('selected',g.dataset.guide===id));select.value=id;
 detail.innerHTML=`<h2>${r.code} · ${esc(r.title)}</h2><div class="english">${esc(r.label)}</div>`+[['input','输入'],['operation','具体操作'],['output','输出'],['next','后续去向'],['example','例子'],['caution','容易误解的地方']].filter(([k])=>r[k]).map(([k,label])=>`<section class="${k==='caution'?'warning':''}"><label>${label}</label><p>${esc(r[k])}</p></section>`).join('')+`<label>实现追溯</label><p class="source">${esc(r.source)}<br>SVG ID: ${esc(r.id)}</p>`;
 const targets=rows.filter(t=>t.id!==id&&r.next.includes(t.code+'「'));if(targets.length){const route=document.createElement('div');route.className='route';targets.forEach(t=>{const b=document.createElement('button');b.textContent='查看 '+t.code;b.onclick=()=>show(t.id);route.appendChild(b)});detail.appendChild(route)}
}
document.querySelectorAll('g[data-guide]').forEach(g=>{g.addEventListener('click',()=>show(g.dataset.guide));g.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();show(g.dataset.guide)}})});
select.onchange=()=>show(select.value);document.getElementById('search').oninput=e=>populate(e.target.value);
document.getElementById('previous').onclick=()=>show(rows[(rows.findIndex(r=>r.id===active)+rows.length-1)%rows.length].id);document.getElementById('next').onclick=()=>show(rows[(rows.findIndex(r=>r.id===active)+1)%rows.length].id);
const svg=document.querySelector('.canvas svg');let width=1800;function resize(w){width=Math.max(600,Math.min(6000,w));svg.style.width=width+'px';}
document.getElementById('smaller').onclick=()=>resize(width/1.2);document.getElementById('larger').onclick=()=>resize(width*1.2);document.getElementById('fit').onclick=()=>resize(document.getElementById('canvas').clientWidth-36);
document.getElementById('locate').onclick=()=>document.getElementById(active).scrollIntoView({block:'center',inline:'center',behavior:'smooth'});
populate();show(rows[0].id);window.refereeGuide={show,populate,getActive:()=>active,rows};
</script></html>'''.replace('__SVG__',svg_inline).replace('__DATA__',encoded)
page=page.replace('href="index.html"','href="'+SVG.name+'"')
(R/'workflow_explained.html').write_text(page)
print('Wrote explanations for',len(rows),'controls:',dict(counts))
