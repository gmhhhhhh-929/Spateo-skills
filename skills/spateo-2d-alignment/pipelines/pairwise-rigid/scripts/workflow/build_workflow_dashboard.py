#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import json
from pathlib import Path


STATE_ORDER = ["slice_level_baseline", "current_review", "current_accepted", "final"]


def load_records(project_root: Path) -> list[tuple[str, dict]]:
    records: list[tuple[str, dict]] = []
    for path in sorted((project_root / "steps").glob("[0-9][0-9][0-9]_*/record.json")):
        try:
            obj = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        records.append((str(path.relative_to(project_root)), obj))
    return records


def load_states(project_root: Path) -> list[tuple[str, dict]]:
    out = []
    for name in STATE_ORDER:
        path = project_root / "states" / f"{name}.yaml"
        if not path.exists():
            continue
        try:
            out.append((f"states/{name}.yaml", json.loads(path.read_text(encoding="utf-8"))))
        except Exception:
            continue
    return out


def link(label: str, href: str) -> str:
    return f'<a href="{html.escape(href)}">{html.escape(label)}</a>'


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a compact chronological spatial workflow dashboard.")
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-html", required=True, type=Path)
    parser.add_argument("--title", default="Spatial alignment workflow")
    args = parser.parse_args()

    state_rows = []
    for state_rel, state in load_states(args.project_root):
        artifacts = state.get("artifacts") or {}
        state_links = [link("state", state_rel)]
        for key in ["clean_coordinates", "viewer", "recipe_chain", "validation", "manifest", "audit_points"]:
            item = artifacts.get(key) or {}
            href = item.get("path") if isinstance(item, dict) else ""
            if href:
                state_links.append(link(key, str(href)))
        clean = (artifacts.get("clean_coordinates") or {}).get("path", "") if isinstance(artifacts.get("clean_coordinates"), dict) else ""
        state_rows.append(
            "<tr>"
            f"<td><code>{html.escape(str(state.get('state', '')))}</code></td>"
            f"<td>{html.escape(str(state.get('status', '')))}</td>"
            f"<td>{html.escape(str(state.get('step', '')))}</td>"
            f"<td>{html.escape(str(state.get('summary', '')))}</td>"
            f"<td>{html.escape(str(state.get('row_count') or ''))}</td>"
            f"<td><code>{html.escape(str(clean))}</code></td>"
            f"<td>{' '.join(state_links)}</td>"
            "</tr>"
        )

    rows = []
    for i, (record_rel, record) in enumerate(load_records(args.project_root), start=1):
        step = Path(record_rel).parent.as_posix()
        outputs = record.get("outputs") or {}
        artifacts = record.get("artifacts") or []
        links = [link("record", record_rel)]
        for key, value in outputs.items():
            if str(value).endswith((".html", ".json", ".csv", ".yaml", ".yml")):
                links.append(link(str(key), f"{step}/{value}"))
        for artifact in artifacts:
            if str(artifact).endswith((".html", ".json", ".csv", ".yaml", ".yml")):
                links.append(link(Path(str(artifact)).name, f"{step}/{artifact}"))
        rows.append(
            "<tr>"
            f"<td>{i}</td>"
            f"<td><code>{html.escape(step)}</code></td>"
            f"<td>{html.escape(str(record.get('summary', '')))}</td>"
            f"<td>{html.escape(str(record.get('status', '')))}</td>"
            f"<td>{html.escape(str(record.get('validator_status', '')))}</td>"
            f"<td>{' '.join(links)}</td>"
            "</tr>"
        )
    doc = "\n".join(
        [
            "<!doctype html><html><head><meta charset='utf-8'>",
            f"<title>{html.escape(args.title)}</title>",
            "<style>body{font:14px -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;margin:24px;color:#17202a}table{border-collapse:collapse;width:100%;max-width:1500px}td,th{border-bottom:1px solid #d8dee6;padding:8px 10px;text-align:left;vertical-align:top}th{background:#f5f7fa}code{font-size:12px;background:#f4f6f8;padding:1px 4px;border-radius:4px}a{color:#0757a6;text-decoration:none}</style>",
            "</head><body>",
            f"<h1>{html.escape(args.title)}</h1>",
            "<h2>Current States</h2>",
            "<table><thead><tr><th>state</th><th>status</th><th>step</th><th>summary</th><th>rows</th><th>clean coordinates</th><th>links</th></tr></thead><tbody>",
            *(state_rows or ["<tr><td colspan='7'>No states found.</td></tr>"]),
            "</tbody></table>",
            "<h2>Operations</h2>",
            "<table><thead><tr><th>#</th><th>step</th><th>action</th><th>status</th><th>validator</th><th>links</th></tr></thead><tbody>",
            *rows,
            "</tbody></table></body></html>",
        ]
    )
    args.output_html.parent.mkdir(parents=True, exist_ok=True)
    args.output_html.write_text(doc, encoding="utf-8")
    print(str(args.output_html))


if __name__ == "__main__":
    main()
