#!/usr/bin/env python3
"""
faln2html.py – render a FASTA alignment as a standalone HTML viewer.

Supported javascript based MSA visualization tools:
  • reactmsa
  • proseqviewer

Examples:
  python faln2html.py -i aln.fa -o out.html --tool reactmsa
  cat aln.fa | python aln2html.py -i - -o out.html --tool reactmsa --no-consensus
"""

import argparse
import json
import sys
import gzip
from pathlib import Path

# ---------------- FASTA Validation ----------------

def validate_fasta_string(fasta_text):
    """Simple FASTA validator: ensures proper headers and sequence lines."""
    lines = fasta_text.splitlines()
    lines = [l for l in lines if l.strip()]
    if not lines:
        raise ValueError("FASTA is empty.")
    if not lines[0].startswith(">"):
        raise ValueError("First line must start with '>'.")
    current_header = None
    has_seq = False
    for i, line in enumerate(lines, 1):
        if line.startswith(">"):
            if current_header and not has_seq:
                raise ValueError("Header at line %d has no sequence lines." % (i - 1))
            current_header = line
            has_seq = False
        else:
            if not all(c.isalpha() or c in "-.*?" for c in line.strip()):
                raise ValueError("Invalid character in line %d: %s" % (i, line))
            has_seq = True
    if not has_seq:
        raise ValueError("Final header has no sequence lines.")



# ---------------- HTML Writer: react-msa ----------------

def write_reactmsa_html(fasta_text, out_html):
    """Write standalone react-msa HTML file."""

    fasta_js = json.dumps(fasta_text)

    html = f"""
<html>
  <head>
    <script
      crossorigin
      src="https://unpkg.com/react-msaview/bundle/index.js"
    ></script>
  </head>
  <body>
    <div id="root" />

    <script>
      const {{ React, createRoot, MSAView, MSAModelF }} = window.ReactMSAView


const fasta = {fasta_js};

const model = MSAModelF().create({{
        type: 'MsaView',
        data: {{ msa: fasta }},
      }})

      model.setColorSchemeName('percent_identity_dynamic')
      model.setBgColor(true)

      // choose MSA width, calculate width of div/rendering area if needed beforehand
      model.setWidth(1300)
      // model.setHeight(800)   // sets vertical extent (depends on implementation)

      const root = createRoot(document.getElementById('root'))
      root.render(React.createElement(MSAView, {{ model }}))
    </script>
  </body>
</html>

    """
    out_path = Path(out_html).resolve()
    out_path.write_text(html, encoding="utf-8")
    return str(out_path)


##########################################################




# ---------------- HTML Writer: ProSeqViewer ----------------

def make_proseqviewer_html(fasta_str, output_path="alignment_view.html", title="ProSeqViewer — FASTA Parsing Minimal", show_consensus=True):
    """Write standalone ProSeqViewer HTML file."""
    fasta_clean = fasta_str.strip().replace("`", "\\`")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <link rel="stylesheet" type="text/css"
        href="https://rawgithub.com/BioComputingUP/ProSeqViewer/master/dist/assets/proseqviewer.css">
  <script src="https://rawgithub.com/BioComputingUP/ProSeqViewer/master/dist/sqv-bundle.js"></script>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 0; }}
    h1 {{ text-align:center; font-size: 18px; margin: 16px 0; }}
    #psv {{
      max-width: 1500px;
      margin: 0 auto 50px;
      white-space: nowrap;
      overflow-x: auto;
      overflow-y: hidden;
      padding-bottom: 16px;
    }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <div id="psv"></div>
  <script>
    const fastaAlignment = `
{fasta_clean}
`.trim();

    function parseFastaAlignment(text) {{
      const records = [];
      let label = null, seq = [];
      const flush = () => {{
        if (label) {{
          const s = seq.join('').replace(/\\s+/g, '');
          records.push({{ label, sequence: s }});
        }}
        label = null; seq = [];
      }};
      for (const rawLine of text.split(/\\r?\\n/)) {{
        const line = rawLine.trim();
        if (!line) continue;
        if (line.startsWith('>')) {{ flush(); label = line.slice(1).trim(); }}
        else seq.push(line.toUpperCase());
      }}
      flush();
      if (records.length === 0) throw new Error('No FASTA records found.');
      const maxLen = Math.max(...records.map(r => r.sequence.length));
      const minLen = Math.min(...records.map(r => r.sequence.length));
      if (minLen !== maxLen) {{
        for (const r of records) {{
          if (r.sequence.length < maxLen) r.sequence += '.'.repeat(maxLen - r.sequence.length);
        }}
      }}
      return records.map((r, i) => ({{ id: i + 1, label: r.label, sequence: r.sequence }}));
    }}

    let sequences = [];
    try {{
      sequences = parseFastaAlignment(fastaAlignment);
    }} catch (e) {{
      console.error(e);
      const errDiv = document.createElement('div');
      errDiv.textContent = 'FASTA parse error: ' + e.message;
      errDiv.style.color = 'crimson';
      errDiv.style.textAlign = 'center';
      document.body.prepend(errDiv);
    }}

    const nucleotideColors = {{A:'#64B964',C:'#4C9BE8',G:'#F0C93D',T:'#E15759',U:'#E15759','N':'#AAAAAA','-':'#FFFFFF','.':'#FFFFFF'}};
    const consensus = {{ dotThreshold: 90, color: nucleotideColors }};
    const options = {{
      chunkSize: 50,
      sequenceColor: nucleotideColors,
      lateralIndexes: false,
      upperIndexes: false,
      showConsensus: {str(show_consensus).lower()},
      lineHeight: 1.2,
      fontSize: '16px'
    }};
    const psv = new ProSeqViewer('psv');
    psv.draw({{ sequences, options, consensus, icons: [] }});
  </script>
</body>
</html>"""

    out_path = Path(output_path).expanduser().resolve()
    out_path.write_text(html, encoding="utf-8")
    return str(out_path)

##########################################################



# ---------------- HTML Writer: MSABrowser ----------------
# Not fully working, large alignments do not load fully...
def write_msabrowser_html(fasta_text, out_html, title="MSABrowser | Inline FASTA",
                          has_consensus=True, pin_version="v1.1", scale=0.5):
    """Write standalone MSABrowser HTML file."""
    if not fasta_text.endswith("\n"):
        fasta_text += "\n"

    fasta_js = json.dumps(fasta_text)
    css_href = f"https://cdn.jsdelivr.net/gh/msabrowser/msabrowser@{pin_version}/css/style.css"
    js_href = f"https://cdn.jsdelivr.net/gh/msabrowser/msabrowser@{pin_version}/javascript/msabrowser.js"
    inv = 100 / max(0.0001, scale)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <title>{title}</title>
  <link rel="stylesheet" href="{css_href}" />
  <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
  <script src="https://html2canvas.hertzen.com/dist/html2canvas.min.js"></script>
  <script src="{js_href}"></script>
  <style>
    body {{ font-family: system-ui, sans-serif; margin: 16px; }}
    #msaZoomWrap {{ transform-origin: top left; transform: scale({scale}); width: {inv:.2f}%; }}
    #MSABrowserDemo {{ margin-top: 12px; }}
  </style>
</head>
<body>
  <h1>{title}</h1>
  <div id="msaZoomWrap"><section id="MSABrowserDemo"></section></div>
  <script>
    const FASTA = {fasta_js};

    const annotations = [];
    const alterations = [];

    $(function () {{
      const viewer = new MSABrowser({{
        id: "MSABrowserDemo",
        msa: MSAProcessor({{
          fasta: FASTA,
          hasConsensus: {str(has_consensus).lower()}
        }}),
        annotations: annotations,
        alterations: alterations,
        colorSchema: "nucleotide"
      }});
    }});
  </script>
</body>
</html>"""
    out_path = Path(out_html).resolve()
    out_path.write_text(html, encoding="utf-8")
    return str(out_path)

##########################################################



# ---------------- Input Handling ----------------

def read_text_auto(path_or_dash):
    """Read from file, .gz, or stdin."""
    if path_or_dash == "-" or (path_or_dash is None and not sys.stdin.isatty()):
        data = sys.stdin.buffer.read()
        try:
            return data.decode("utf-8")
        except UnicodeDecodeError:
            return data.decode("latin-1")
    if not path_or_dash:
        sys.exit("Error: no input provided. Use -i FILE, -i -, or pipe FASTA on stdin.")
    p = Path(path_or_dash)
    if not p.exists():
        sys.exit(f"Error: file not found: {p}")
    if p.suffix == ".gz":
        with gzip.open(p, "rt", encoding="utf-8", errors="replace") as f:
            return f.read()
    return p.read_text(encoding="utf-8", errors="replace")


# ---------------- CLI ----------------

def main():
    parser = argparse.ArgumentParser(description="Embedd a FASTA alignment in a HTML file using a javascript based viewer")
    parser.add_argument("-i", "--input", help="FASTA file path or '-' for STDIN.")
    parser.add_argument("-o", "--out", required=True, help="Output HTML file path.")
    parser.add_argument("--tool",
                        choices=["reactmsa", "proseqviewer"], # "msabrowser",
                        default="reactmsa",
                        help="Which viewer to use (default: msabrowser).")
    parser.add_argument("--no-consensus", action="store_true",
                        help="Disable consensus row (applies to proseqviewer).")
    #parser.add_argument("--title", default="Alignment Viewer", help="HTML page title.")

                        # help="Disable consensus row (applies to msabrowser and proseqviewer).")
    # parser.add_argument("--scale", type=float, default=0.5,
    #                     help="Zoom scale (MSABrowser only).")

    args = parser.parse_args()
    fasta_text = read_text_auto(args.input)

    try:
        validate_fasta_string(fasta_text)
    except ValueError as e:
        sys.exit(f"FASTA validation error: {e}")

    has_consensus = not args.no_consensus

    if args.tool == "reactmsa":
        out = write_reactmsa_html(fasta_text, args.out)

    elif args.tool == "proseqviewer":
        out = make_proseqviewer_html(
            fasta_text, args.out, title=args.title,
            show_consensus=has_consensus)

    elif args.tool == "msabrowser":
        out = write_msabrowser_html(
            fasta_text, args.out, title=args.title,
            has_consensus=has_consensus, scale=args.scale)

    else:
        sys.exit("Unknown tool option.")

    print(f"Wrote HTML file to: {out}", file=sys.stderr)

if __name__ == "__main__":
    main()
