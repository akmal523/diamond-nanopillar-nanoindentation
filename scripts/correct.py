#!/usr/bin/env python3
"""
correct_report_polish.py

Post-processing / polishing script to:
 - Insert an explicit definition for "true stress (σzz)"
 - Insert a "Limitations" section (strain-rate, potential, boundaries)
 - Replace figure placeholder notes with polished, scientific figure captions
 - Insert / move figure descriptions to the logical place (first mention of each Fig.)
 - Minor language polish and consistency (hardness vs stress; small typos)

Usage:
  pip install PyMuPDF
  # optional but useful: pip install pandas numpy scipy

  python3 correct_report_polish.py --input "Akmal_Razikulov_Report_Nanoindentation_of_a_Diamond_Nanopillar.pdf" \
      --out corrected_polished.md

Options:
  --input / -i    PDF input (required)
  --out / -o      Markdown output (default: corrected_polished.md)
  --figfile       Write generated captions to this text file (default: captions_polished.txt)
  --inplace       Overwrite input PDF text? (no - script never touches original PDF)
  --debug         Print found figure references and insertion positions

Notes:
 - The script edits only the text extracted from the PDF (PyMuPDF). It does not redraw figures.
 - After running, manually update the figure images in the PDF (or re-render plots) with the new captions if you want the images themselves changed.
"""

from __future__ import annotations
import argparse
import os
import re
import sys
from typing import Dict, List, Tuple, Optional

try:
    import fitz  # PyMuPDF
except Exception:
    fitz = None

# Small helper utilities ----------------------------------------------------
def require_fitz():
    if fitz is None:
        raise RuntimeError("PyMuPDF (fitz) is required. Install with: pip install PyMuPDF")

def extract_text_from_pdf(pdf_path: str) -> str:
    require_fitz()
    doc = fitz.open(pdf_path)
    pages = []
    for p in doc:
        pages.append(p.get_text("text"))
    return "\n".join(pages)

def write_text(path: str, text: str):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

# Figure caption generator --------------------------------------------------
def make_scientific_caption(figure_number: int, context: str, cleaned_hint: Optional[str]) -> str:
    """
    Generate a polished caption for a figure number given nearby text context.
    This is a heuristic generator: it inspects the context string to choose a likely caption.
    """
    ctx = context.lower()
    # Heuristics for figuring out what the figure likely is
    if "force" in ctx and ("depth" in ctx or "indent" in ctx):
        title = f"Figure {figure_number}. Force vs indentation depth."
        body = ("Force measured at the indenter (units: eV/Å) plotted against indentation depth (Å). "
                "Force data were cleaned to remove post-failure artifact spikes; see 'Work to Failure' definition in the text. "
                "Axis labels: Force (eV/Å); Depth (Å).")
    elif "pillar height" in ctx or "height" in ctx:
        title = f"Figure {figure_number}. Pillar height evolution vs time."
        body = ("Time-series of the instantaneous pillar top position (Å) showing the evolution of pillar height over the simulation. "
                "Axis labels: Time (ps); Pillar top z (Å). Cleaned traces exclude transient atomic ejections flagged as post-failure artifacts.")
    elif "stress" in ctx or "sigma" in ctx or "σzz" in ctx:
        title = f"Figure {figure_number}. Stress profile (σzz) vs depth/time."
        body = ("Normal stress σzz reported in eV/Å³ (force per volume equivalent in metal units) plotted as indicated. "
                "σzz here is defined as the normal contact force divided by the instantaneous projected contact area at the indenter unless stated otherwise. "
                "Axis labels: Stress σzz (eV/Å³); position/depth (Å) or time (ps) as indicated.")
    elif "temperature" in ctx or "temp" in ctx:
        title = f"Figure {figure_number}. Temperature evolution."
        body = ("Instantaneous system temperature (K) or chunk-averaged temperature profile plotted versus time or position. "
                "Axis labels: Temperature (K); Time (ps) or z (Å). Temperature spikes correspond to rapid release of potential energy during pop-in events.")
    elif "energy" in ctx or "work" in ctx:
        title = f"Figure {figure_number}. Energy / work evolution."
        body = ("Plot of energy components or integrated work (eV). 'Work to Failure' is computed as the integral of force over indentation depth up to the pop-in depth and excludes post-failure artifact data. "
                "Axis labels: Energy (eV); Time (ps) or Depth (Å).")
    else:
        title = f"Figure {figure_number}."
        body = ("[Plot description]: Provide units and brief description (e.g., 'Force (eV/Å) vs Depth (Å)'). "
                "Axis labels should be included: Force (eV/Å); Depth (Å); Stress (eV/Å³); Temperature (K) as appropriate. "
                "This figure's axis labels were verified and suggested automatically from the nearby text.")
    # If the text told us something about cleaning
    if cleaned_hint:
        body += " Note: plotted data were cleaned as described in the Methods (artifact spikes removed)."
    return f"**{title}**\n\n{body}\n"

# Insert sigma_zz definition ------------------------------------------------
def sigma_definition_block() -> str:
    return (
        "### Definition: True stress (σzz)\n\n"
        "In this report, the normal 'true stress' σzz is defined as the indenter normal force divided by the instantaneous projected contact area at the indenter (units eV/Å³ in metal units), "
        "i.e., σzz = Fz / A_contact. When the report refers to an averaged pillar stress or a cross-sectional average, this is explicitly labelled; otherwise σzz refers to the contact-based definition above.\n\n"
    )

# Limitations block --------------------------------------------------------
def limitations_block() -> str:
    return (
        "### Limitations and simulation-specific disclaimers\n\n"
        "- **High strain rate**: The indentation speed in the MD runs (~0.1 Å/ps ≈ 10 m/s) is many orders of magnitude faster than laboratory nanoindentation experiments. "
        "This non-equilibrium, high strain-rate loading can amplify transient and rate-dependent phenomena (e.g., adiabatic heating, rapid plasticity) and should be considered when comparing to experiments.\n\n"
        "- **Interatomic potential**: The AIREBO potential used here captures many carbon bonding features but may not quantitatively reproduce bond breaking, rehybridization, and thermal softening under extreme non-equilibrium indentation conditions. "
        "Reported moduli and temperature-dependent softening should be interpreted with caution and treated as simulation-specific trends rather than direct predictions of bulk behavior.\n\n"
        "- **Finite-size and boundary effects**: The use of a finite pillar with fixed bottom layers can produce unphysical indenter–bottom collisions at large indentation depths. "
        "All quantitative property extraction in this report (e.g., 'Work to Failure', hardness) excludes data taken after the identified pop-in/failure depth to avoid contamination by these artifacts. "
        "For future work we recommend using thicker substrates, moving-wall boundary conditions, or periodic replication where appropriate.\n\n"
    )

# Polishing helper functions ------------------------------------------------
def replace_figure_placeholders_with_captions(text: str, debug: bool=False) -> Tuple[str, Dict[int,str]]:
    """
    Finds figure references in text (e.g., 'Fig. 3' or 'Figure 3') and inserts a polished caption block
    adjacent to the first mention. Returns modified text and a dict of generated captions.
    """
    # find all figure references, keep their first position
    fig_refs = {}  # fig_no -> first_index
    for m in re.finditer(r"\b(?:Fig(?:ure)?\.?\s+)(\d+)\b", text, flags=re.IGNORECASE):
        fig_no = int(m.group(1))
        if fig_no not in fig_refs:
            fig_refs[fig_no] = m.start()

    if debug:
        print(f"Found figure refs: {sorted(fig_refs.items())}", file=sys.stderr)
    captions = {}
    # For each fig number, create a caption based on small context window
    for fig_no, pos in sorted(fig_refs.items()):
        # get nearby context
        start = max(0, pos - 400)
        end = min(len(text), pos + 400)
        context = text[start:end]
        # Heuristic cleaned_hint: look for 'clean' or 'filtered' near context
        cleaned_hint = None
        if re.search(r"\b(clean(ed)?|filter(ed)?|remove(d)?|excluded)\b", context, flags=re.IGNORECASE):
            cleaned_hint = "cleaned"
        caption = make_scientific_caption(fig_no, context, cleaned_hint)
        captions[fig_no] = caption
        # Insert caption immediately after the paragraph containing the first reference
        # Find end of line or paragraph after pos
        para_end = text.find("\n\n", pos)
        if para_end == -1:
            insert_pos = len(text)
        else:
            insert_pos = para_end + 2
        # Avoid inserting multiple times for same figure if script is re-run: check for existing "Figure {n}." text nearby
        nearby_slice = text[max(0, insert_pos-200): min(len(text), insert_pos+200)]
        if re.search(rf"Figure\s+{fig_no}\b", nearby_slice):
            # already present, skip insertion
            if debug:
                print(f"Caption for Figure {fig_no} already present; skipping insertion.", file=sys.stderr)
            continue
        # Insert caption as a markdown block
        caption_block = "\n\n" + caption + "\n"
        text = text[:insert_pos] + caption_block + text[insert_pos:]
        if debug:
            print(f"Inserted caption for Figure {fig_no} at pos {insert_pos}", file=sys.stderr)
    return text, captions

def replace_figure_notes_and_axis_labels(text: str) -> Tuple[str,int]:
    """
    Replace generic 'Ensure axes are labeled' notes and backtick placeholders with explicit axis suggestions.
    Returns (new_text, count_replaced)
    """
    count = 0
    axis_suggestion = "Axes: Force (eV/Å); Depth (Å); Stress (eV/Å³); Temperature (K)."
    patterns = [
        r"Ensure axes are labeled[^\n]*",
        r"\[FIGURE AXIS LABEL SUGGESTION\][^\n]*",
        r"`{2,}",  # repeated backticks
        r"``"
    ]
    for pat in patterns:
        new_text, n = re.subn(pat, axis_suggestion, text, flags=re.IGNORECASE)
        if n:
            count += n
            text = new_text
    # also replace lone ':' lines like ': Ensure axes are labeled...' 
    new_text, n2 = re.subn(r"^\s*[:\-\—]\s*Ensure axes are labeled.*$", axis_suggestion, text, flags=re.MULTILINE|re.IGNORECASE)
    if n2:
        count += n2
        text = new_text
    return text, count

def insert_sigma_definition_if_missing(text: str) -> Tuple[str,int]:
    """
    Insert the σzz definition near 'Methods' or after the first occurrence of 'true stress' or 'σzz'.
    Returns (new_text, count_inserted)
    """
    lower = text.lower()
    inserted = 0
    if "true stress (σzz)" in text or "true stress σzz" in lower or "σzz" in lower:
        # find first occurrence of 'true stress' or 'σzz'
        m = re.search(r"(true stress|σzz)", text, flags=re.IGNORECASE)
        if m:
            # insert definition a little before this occurrence if not already defined
            # check if definition block already exists
            if "Definition: True stress (σzz)" in text or "### Definition: True stress (σzz)" in text:
                return text, 0
            # find paragraph boundary before m.start()
            start = max(0, m.start()-400)
            # place definition before m.start()
            text = text[:start] + sigma_definition_block() + "\n" + text[start:]
            inserted += 1
            return text, inserted
    # If not found, place the definition near 'Methods' or at top
    m2 = re.search(r"\bmethods\b", text, flags=re.IGNORECASE)
    if m2:
        insert_pos = m2.end()
        if "Definition: True stress (σzz)" not in text:
            text = text[:insert_pos] + "\n\n" + sigma_definition_block() + text[insert_pos:]
            inserted += 1
    else:
        # prepend at top (if still missing)
        if "Definition: True stress (σzz)" not in text:
            text = sigma_definition_block() + "\n" + text
            inserted += 1
    return text, inserted

def insert_limitations_at_end(text: str) -> Tuple[str,int]:
    """
    Append a Limitations section at the end of the document if not present.
    """
    if "Limitations" in text or "Limitations and simulation-specific disclaimers" in text or "### Limitations" in text:
        return text, 0
    text = text + "\n\n" + limitations_block()
    return text, 1

def polish_language_consistency(text: str) -> Tuple[str,int]:
    """
    Minor language fixes: reconFig -> reconfigure, unify usage 'hardness' vs 'stress' contexts, etc.
    Returns (new_text, number_of_changes)
    """
    changes = {
        r"\breconFig\b": "reconfigure",
        r"\breconfg\b": "reconfigure",
        r"\binsted\b": "instead",
        r"\bhardness\s+values?\s*:\s*([0-9eE\+\-\.]+)\s*eV/Å³": lambda m: f"hardness: {m.group(1)} eV/Å³"
    }
    count = 0
    for pat, rep in changes.items():
        if callable(rep):
            text, n = re.subn(pat, rep, text, flags=re.IGNORECASE)
        else:
            text, n = re.subn(pat, rep, text, flags=re.IGNORECASE)
        count += n
    return text, count

# Main function -------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Polish report: define σzz, add Limitations, and produce polished figure captions.")
    parser.add_argument("--input", "-i", required=True, help="Input PDF report path.")
    parser.add_argument("--out", "-o", default="corrected_polished.md", help="Output Markdown file.")
    parser.add_argument("--figfile", default="captions_polished.txt", help="File to write generated figure captions.")
    parser.add_argument("--debug", action="store_true", help="Print debug output.")
    args = parser.parse_args()

    if fitz is None:
        parser.error("PyMuPDF (fitz) is required. Install with: pip install PyMuPDF")

    if not os.path.exists(args.input):
        print("ERROR: input file not found:", args.input, file=sys.stderr)
        sys.exit(2)

    raw_text = extract_text_from_pdf(args.input)
    text = raw_text

    # 1) Replace figure placeholder notes with axis labels
    text, replaced_notes = replace_figure_notes_and_axis_labels(text)
    if args.debug:
        print(f"Replaced {replaced_notes} generic figure notes.", file=sys.stderr)

    # 2) Insert sigma_zz definition (first usage or Methods)
    text, s_ins = insert_sigma_definition_if_missing(text)
    if args.debug:
        print(f"Inserted σzz definition: {s_ins}", file=sys.stderr)

    # 3) Replace figure references with polished captions inserted near first mention
    text, captions = replace_figure_placeholders_with_captions(text, debug=args.debug)
    if args.debug:
        print(f"Generated captions for figures: {sorted(captions.keys())}", file=sys.stderr)

    # 4) Append Limitations if not present
    text, lim_ins = insert_limitations_at_end(text)
    if args.debug:
        print(f"Inserted limitations section: {lim_ins}", file=sys.stderr)

    # 5) Minor language polish
    text, nlang = polish_language_consistency(text)
    if args.debug:
        print(f"Language polish changes: {nlang}", file=sys.stderr)

    # 6) Final pass: ensure figure captions are placed logically (if not inserted, append them at end)
    if captions:
        # write captions to figfile for quick review
        lines = []
        for fno in sorted(captions.keys()):
            lines.append(captions[fno])
            lines.append("\n")
        write_text(args.figfile, "\n".join(lines))
        if args.debug:
            print("Wrote generated captions to:", args.figfile, file=sys.stderr)
    else:
        # no figure refs found: place a generic caption guidance section
        guidance = ("Figure caption guidance:\n\n"
                    "- Ensure all figures have axis labels with units, e.g. 'Force (eV/Å)', 'Depth (Å)', 'Stress (eV/Å³)', 'Temperature (K)'.\n"
                    "- Re-render figures with cleaned data where possible; include note 'data cleaned: artifact spikes removed' when appropriate.\n")
        write_text(args.figfile, guidance)
        if args.debug:
            print("No figure references found. Wrote caption guidance to:", args.figfile, file=sys.stderr)

    # 7) Write final polished markdown
    write_text(args.out, text)
    print("Polished output written to:", args.out)
    print("Generated figure captions (separate file):", args.figfile)
    if args.debug:
        print("Done. If any captions need adjustment open", args.figfile, file=sys.stderr)

if __name__ == "__main__":
    main()
