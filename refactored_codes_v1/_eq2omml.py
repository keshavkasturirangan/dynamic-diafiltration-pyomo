#!/usr/bin/env python3
"""LaTeX math -> OOXML OMML (native, editable) via pandoc, for PPTX injection.

Caches by md5 so repeated equations are converted once.
"""
import hashlib
import re
import subprocess
import zipfile
from pathlib import Path

_CACHE = Path("/tmp/_omml_cache")
_CACHE.mkdir(exist_ok=True)

_OMATH_RE = re.compile(r"<m:oMath>.*?</m:oMath>", re.S)


def latex_to_omath(latex: str) -> str:
    """Return the inner <m:oMath>...</m:oMath> string (m: namespace) for a LaTeX
    math expression (no surrounding $)."""
    key = hashlib.md5(latex.encode()).hexdigest()
    cached = _CACHE / f"{key}.xml"
    if cached.exists():
        return cached.read_text()
    md = _CACHE / f"{key}.md"
    docx = _CACHE / f"{key}.docx"
    md.write_text(f"$${latex}$$\n")
    subprocess.run(["pandoc", "-f", "markdown", "-t", "docx", "-o", str(docx), str(md)],
                   check=True, capture_output=True)
    with zipfile.ZipFile(docx) as z:
        xml = z.read("word/document.xml").decode()
    m = _OMATH_RE.search(xml)
    if not m:
        raise RuntimeError(f"no oMath produced for: {latex}\n{xml[:400]}")
    omath = m.group(0)
    cached.write_text(omath)
    return omath


if __name__ == "__main__":
    import sys
    print(latex_to_omath(sys.argv[1] if len(sys.argv) > 1
                         else r"J_w = L_p\left(\Delta P - \sigma\,\Delta\pi\right)")[:800])
