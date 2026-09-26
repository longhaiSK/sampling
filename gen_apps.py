#!/bin/sh
""":"
# 1. Try Unix venv
PYTHON_EXEC="$(dirname "$0")/.venv/bin/python"

# 2. If not found, try Windows venv
if [ ! -f "$PYTHON_EXEC" ]; then
    PYTHON_EXEC="$(dirname "$0")/.venv/Scripts/python"
fi

# 3. If still not found, try system PATH
if [ ! -f "$PYTHON_EXEC" ]; then
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_EXEC="python3"
    elif command -v python >/dev/null 2>&1; then
        PYTHON_EXEC="python"
    else
        echo "gen_gallery: Neither virtualenv nor system Python found." >&2
        exit 1
    fi
fi

exec "$PYTHON_EXEC" "$0" "$@"
"""

# Regenerates gallery.qmd (the "Figures, Tables, and Apps" index) straight from

# Regenerates the standalone Shinylive app pages and their index,
# shinyliveapps_sampling.qmd, straight from the app blocks already in the
# chapter .qmd files, so the archive can never drift out of sync with the book.
# Run it after adding, removing or editing any shinylive app:
#
#   ./gen_apps.py && quarto render shinyliveapps_sampling.qmd app-*.qmd
#
# Why separate pages rather than one page of tabs: every app boots its own webR
# runtime (tens of MB), so six on one page would be unusable. One app per page
# keeps each load to a single runtime; the top navigation makes switching cheap.
#
# These pages are invisible to the book build -- a Quarto *book* project renders
# only the files listed in _quarto.yml's chapters/appendices -- so they carry
# their own theme and full-width layout without touching the book at all.

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
QUARTO_YML = ROOT / "_quarto.yml"
INDEX = ROOT / "shinyliveapps_sampling.qmd"
CSS = ROOT / "shinyliveapps.css"

# slug (page name) and short navigation label for each known app
APP_META = {
    "fig-app-srs":       ("srs",        "SRS"),
    "fig-strat-app":     ("stratified", "Stratified"),
    "fig-ratio-app":     ("ratio",      "Ratio &amp; regression"),
    "fig-poststrat-app": ("poststrat",  "Post-stratification"),
    "fig-cluster-app":   ("cluster",    "Cluster"),
    "fig-ups-app":       ("ups",        "UPS"),
}

APP_DIV_RE = re.compile(r'^::: \{#(fig-[\w-]+)\}[ \t]*\n(.*?)^:::[ \t]*$', re.S | re.M)
CODE_RE    = re.compile(r'^```\{shinylive-r\}[ \t]*\n(.*?)^```[ \t]*$', re.S | re.M)
CAPTION_RE = re.compile(r'^\*\*(.+?)\*\*[ \t]*$', re.M)
TITLE_RE   = re.compile(r'^#\s+(.+)$', re.M)


def chapter_order():
    text = QUARTO_YML.read_text(encoding="utf-8")
    out, in_list = [], False
    for line in text.splitlines():
        if re.match(r'^\s*(chapters|appendices):\s*$', line):
            in_list = True
            continue
        if in_list:
            m = re.match(r'^\s*-\s*(\S+\.qmd)\s*$', line)
            if m:
                out.append(m.group(1))
                continue
            if line.strip() and not line.startswith((' ', '-', '\t')):
                in_list = False
    return out


def chapter_title(path):
    m = TITLE_RE.search(path.read_text(encoding="utf-8"))
    return m.group(1).strip() if m else path.stem


def find_apps(path):
    """Every fig- div in `path` that contains a shinylive-r block."""
    text = path.read_text(encoding="utf-8")
    found = []
    for m in APP_DIV_RE.finditer(text):
        div_id, body = m.group(1), m.group(2)
        code = CODE_RE.search(body)
        if not code:
            continue
        caps = CAPTION_RE.findall(body)
        caption = caps[-1].strip() if caps else div_id
        found.append({"id": div_id, "code": code.group(1).rstrip("\n"),
                      "caption": caption, "chapter": path.name,
                      "chapter_title": chapter_title(path)})
    return found


def nav_html(apps, active_slug):
    items = ['<a href="shinyliveapps_sampling.html"%s>All apps</a>'
             % (' class="active"' if active_slug is None else "")]
    for a in apps:
        cls = ' class="active"' if a["slug"] == active_slug else ""
        items.append('<a href="app-%s.html"%s>%s</a>' % (a["slug"], cls, a["nav"]))
    items.append('<a class="book" href="index.html">&#8592; back to the book</a>')
    return ('```{=html}\n<nav class="appnav">\n  '
            + "\n  ".join(items) + "\n</nav>\n```\n")


def app_page(app, apps):
    return f"""---
title: "{app['caption']}"
format:
  html:
    theme: flatly
    page-layout: full
    toc: false
    css: shinyliveapps.css
engine: knitr
filters:
  - shinylive
---

{nav_html(apps, app['slug'])}
::: {{.appsource}}
From [{app['chapter_title']}]({Path(app['chapter']).with_suffix('.html')}#{app['id']})
in *Elements of Sampling Survey*.
:::

```{{shinylive-r}}
{app['code']}
```
"""


def index_page(apps):
    cards = []
    for a in apps:
        cards.append(
            f'<a class="appcard" href="app-{a["slug"]}.html">\n'
            f'  <span class="appcard-name">{a["nav"]}</span>\n'
            f'  <span class="appcard-desc">{a["caption"]}</span>\n'
            f'  <span class="appcard-from">{a["chapter_title"]}</span>\n'
            f'</a>'
        )
    return f"""---
title: "Shinylive Apps for Sampling Survey"
format:
  html:
    theme: flatly
    page-layout: full
    toc: false
    css: shinyliveapps.css
engine: markdown
---

{nav_html(apps, None)}
The simulation apps from [Elements of Sampling Survey](index.html), each on its
own page so that only one webR runtime has to start. They run entirely in the
browser --- nothing is sent to a server, and the first load of a page takes a
few seconds while R itself is fetched.

```{{=html}}
<div class="appgrid">
{chr(10).join(cards)}
</div>
```
"""


CSS_TEXT = """/* Shinylive app archive --- deliberately not the book's look.
   Wide layout so an app never has to squeeze itself into a text column. */

main.content, #quarto-content {
  max-width: none !important;
  padding-left: 1.5rem;
  padding-right: 1.5rem;
}
.page-columns { display: block !important; }

/* top navigation */
nav.appnav {
  display: flex;
  flex-wrap: wrap;
  gap: 0.15rem 0.4rem;
  align-items: center;
  margin: 0 0 1.1rem 0;
  padding: 0.45rem 0.2rem;
  border-bottom: 2px solid #18bc9c;
}
nav.appnav a {
  padding: 0.3rem 0.7rem;
  border-radius: 4px;
  text-decoration: none;
  font-weight: 600;
  font-size: 0.92rem;
  color: #2c3e50;
}
nav.appnav a:hover   { background: #e8f6f3; }
nav.appnav a.active  { background: #18bc9c; color: #fff; }
nav.appnav a.book    { margin-left: auto; font-weight: 400; color: #7b8a8b; }

/* provenance line under the title */
.appsource {
  margin: -0.4rem 0 1rem 0;
  font-size: 0.9rem;
  color: #7b8a8b;
}

/* index cards */
.appgrid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
  gap: 1rem;
  margin-top: 1.2rem;
}
a.appcard {
  display: block;
  padding: 0.9rem 1rem;
  border: 1px solid #dfe6e9;
  border-radius: 6px;
  text-decoration: none;
  color: #2c3e50;
  background: #fff;
  transition: border-color .15s, box-shadow .15s;
}
a.appcard:hover { border-color: #18bc9c; box-shadow: 0 2px 8px rgba(0,0,0,.07); }
.appcard-name { display: block; font-weight: 700; font-size: 1.05rem; }
.appcard-desc { display: block; margin-top: .35rem; font-size: .92rem; color: #56666f; }
.appcard-from { display: block; margin-top: .5rem; font-size: .8rem; color: #95a5a6; }
"""


def main():
    apps = []
    for name in chapter_order():
        path = ROOT / name
        if not path.exists():
            continue
        for app in find_apps(path):
            slug, nav = APP_META.get(
                app["id"], (app["id"].replace("fig-", "").replace("-app", ""),
                            app["id"].replace("fig-", "").replace("-app", "")))
            app["slug"], app["nav"] = slug, nav
            apps.append(app)

    if not apps:
        print("gen_apps: no shinylive apps found", file=sys.stderr)
        return 1

    for a in apps:
        (ROOT / f"app-{a['slug']}.qmd").write_text(app_page(a, apps), encoding="utf-8")
    INDEX.write_text(index_page(apps), encoding="utf-8")
    CSS.write_text(CSS_TEXT, encoding="utf-8")

    print(f"gen_apps: wrote {INDEX.name}, {len(apps)} app pages and {CSS.name}")
    for a in apps:
        print(f"   app-{a['slug']}.qmd  <-  {a['chapter']} #{a['id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
