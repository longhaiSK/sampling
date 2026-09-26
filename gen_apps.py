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
INDEX = ROOT / "shinyliveapps_sampling.qmd"
CSS = ROOT / "shinyliveapps.css"

# The app pages ARE the source: each app-<slug>.qmd holds one shinylive app and
# its documentation. This script only keeps the shared furniture in sync -- the
# top navigation on every page, the index page, and the stylesheet -- so adding
# an app means writing app-<slug>.qmd and adding a row to ORDER below.
ORDER = [
    ("srs",        "SRS"),
    ("stratified", "Stratified"),
    ("ratio",      "Ratio &amp; regression"),
    ("poststrat",  "Post-stratification"),
    ("cluster",    "Cluster"),
    ("ups",        "UPS"),
]

FM_RE  = re.compile(r'^---\n(.*?)\n---\n', re.S)
NAV_RE = re.compile(r'<!-- nav:start -->.*?<!-- nav:end -->', re.S)


def front_matter(text, key):
    m = FM_RE.match(text)
    if not m:
        return ""
    v = re.search(rf'^{key}:\s*"?(.*?)"?\s*$', m.group(1), re.M)
    return v.group(1).strip() if v else ""


def nav_html(apps, active):
    items = ['<a href="shinyliveapps_sampling.html"%s>All apps</a>'
             % (' class="active"' if active is None else "")]
    for slug, label in apps:
        items.append('<a href="app-%s.html"%s>%s</a>'
                     % (slug, ' class="active"' if slug == active else "", label))
    items.append('<a class="book" href="index.html">&#8592; back to the book</a>')
    return ("<!-- nav:start -->\n```{=html}\n<nav class=\"appnav\">\n  "
            + "\n  ".join(items) + "\n</nav>\n```\n<!-- nav:end -->")


def index_page(apps, meta):
    cards = "\n".join(
        f'<a class="appcard" href="app-{slug}.html">\n'
        f'  <span class="appcard-name">{label}</span>\n'
        f'  <span class="appcard-desc">{meta[slug]["desc"]}</span>\n'
        f'</a>' for slug, label in apps)
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
few seconds while R itself is fetched. Every page carries the app's
documentation underneath it.

```{{=html}}
<div class="appgrid">
{cards}
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
    apps, meta, missing = [], {}, []
    for slug, label in ORDER:
        page = ROOT / f"app-{slug}.qmd"
        if not page.exists():
            missing.append(page.name)
            continue
        text = page.read_text(encoding="utf-8")
        meta[slug] = {"desc": front_matter(text, "description") or
                              front_matter(text, "title") or slug}
        apps.append((slug, label))

    if missing:
        print("gen_apps: missing pages: " + ", ".join(missing), file=sys.stderr)
    if not apps:
        return 1

    for slug, _ in apps:
        page = ROOT / f"app-{slug}.qmd"
        text = page.read_text(encoding="utf-8")
        nav = nav_html(apps, slug)
        if NAV_RE.search(text):
            text = NAV_RE.sub(lambda _: nav, text, count=1)
        else:
            text = FM_RE.sub(lambda m: m.group(0) + "\n" + nav + "\n", text, count=1)
        page.write_text(text, encoding="utf-8")

    INDEX.write_text(index_page(apps, meta), encoding="utf-8")
    CSS.write_text(CSS_TEXT, encoding="utf-8")
    print(f"gen_apps: refreshed navigation on {len(apps)} app pages, "
          f"plus {INDEX.name} and {CSS.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
