"""
Autonomous Lab - Citation Utilities

Provides DOI-to-BibTeX lookup, title-based paper search, and BibTeX
validation using the CrossRef and PubMed APIs (free, no API key needed).

All functions are designed to be called by the Trainee during paper writing
or by the PI during citation QA.
"""

import json
import re
import urllib.parse
import urllib.request
from typing import Optional


# ---------------------------------------------------------------------------
# CrossRef API helpers
# ---------------------------------------------------------------------------

_CROSSREF_API = "https://api.crossref.org"
_HEADERS = {
    "User-Agent": "AutonomousLab/0.5 (mailto:autolab@example.com)",
    "Accept": "application/json",
}


def _crossref_get(path: str, params: dict | None = None, timeout: int = 15) -> dict:
    """Make a GET request to the CrossRef API."""
    url = f"{_CROSSREF_API}{path}"
    if params:
        url += "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers=_HEADERS)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _format_authors_bibtex(authors: list[dict]) -> str:
    """Format CrossRef author list into BibTeX author string."""
    parts = []
    for a in authors:
        family = a.get("family", "")
        given = a.get("given", "")
        if family and given:
            parts.append(f"{family}, {given}")
        elif family:
            parts.append(family)
    return " and ".join(parts)


def _make_citation_key(item: dict) -> str:
    """Generate a citation key like 'Smith2024keyword'."""
    authors = item.get("author", [])
    first_author = authors[0].get("family", "Unknown") if authors else "Unknown"
    # Clean the name
    first_author = re.sub(r"[^a-zA-Z]", "", first_author)

    year = ""
    date_parts = item.get("published-print", item.get("published-online", {}))
    if date_parts and "date-parts" in date_parts:
        year = str(date_parts["date-parts"][0][0])

    title = item.get("title", [""])[0]
    # Extract first meaningful word from title
    title_words = re.findall(r"[A-Za-z]{4,}", title)
    keyword = title_words[0].lower() if title_words else ""

    return f"{first_author}{year}{keyword}"


def _item_to_bibtex(item: dict) -> str:
    """Convert a CrossRef work item to a BibTeX entry string."""
    entry_type = "article"
    container = item.get("container-title", [""])[0]
    if item.get("type") == "proceedings-article":
        entry_type = "inproceedings"
    elif item.get("type") == "book":
        entry_type = "book"

    key = _make_citation_key(item)
    authors = _format_authors_bibtex(item.get("author", []))
    title = item.get("title", [""])[0]
    year = ""
    date_parts = item.get("published-print", item.get("published-online", {}))
    if date_parts and "date-parts" in date_parts:
        year = str(date_parts["date-parts"][0][0])

    doi = item.get("DOI", "")
    volume = item.get("volume", "")
    issue = item.get("issue", "")
    pages = item.get("page", "")
    publisher = item.get("publisher", "")

    # Protect capitalization in title
    protected_title = re.sub(
        r"\b([A-Z][A-Z]+)\b", lambda m: "{" + m.group(1) + "}", title
    )

    lines = [f"@{entry_type}{{{key},"]
    lines.append(f"  author    = {{{authors}}},")
    lines.append(f"  title     = {{{protected_title}}},")
    if container:
        field = "booktitle" if entry_type == "inproceedings" else "journal"
        lines.append(f"  {field:<9} = {{{container}}},")
    lines.append(f"  year      = {{{year}}},")
    if volume:
        lines.append(f"  volume    = {{{volume}}},")
    if issue:
        lines.append(f"  number    = {{{issue}}},")
    if pages:
        lines.append(f"  pages     = {{{pages}}},")
    if publisher and entry_type == "book":
        lines.append(f"  publisher = {{{publisher}}},")
    if doi:
        lines.append(f"  doi       = {{{doi}}},")
    lines.append("}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def doi_to_bibtex(doi: str) -> str:
    """Look up a DOI via CrossRef and return a BibTeX entry.

    Args:
        doi: A DOI string (e.g., "10.1038/s41586-021-03819-2")

    Returns:
        A formatted BibTeX entry string.

    Raises:
        ValueError: If the DOI cannot be resolved.
    """
    doi = doi.strip().removeprefix("https://doi.org/").removeprefix("http://doi.org/")
    try:
        data = _crossref_get(f"/works/{urllib.parse.quote(doi, safe='')}")
        item = data["message"]
        return _item_to_bibtex(item)
    except Exception as e:
        raise ValueError(f"Could not resolve DOI '{doi}': {e}") from e


def search_papers(
    query: str,
    rows: int = 5,
    sort: str = "relevance",
    filter_from_year: int | None = None,
) -> list[dict]:
    """Search CrossRef for papers matching a query.

    Args:
        query: Search string (title, keywords, etc.)
        rows: Number of results to return (max 20)
        sort: Sort order — "relevance" or "published"
        filter_from_year: Only return papers from this year onward

    Returns:
        List of dicts with keys: title, authors, year, doi, journal, bibtex
    """
    params: dict = {
        "query": query,
        "rows": min(rows, 20),
        "sort": sort,
        "order": "desc",
    }
    if filter_from_year:
        params["filter"] = f"from-pub-date:{filter_from_year}"

    try:
        data = _crossref_get("/works", params)
    except Exception as e:
        return [{"error": f"Search failed: {e}"}]

    results = []
    for item in data.get("message", {}).get("items", []):
        title = item.get("title", [""])[0]
        authors = _format_authors_bibtex(item.get("author", []))
        year = ""
        date_parts = item.get("published-print", item.get("published-online", {}))
        if date_parts and "date-parts" in date_parts:
            year = str(date_parts["date-parts"][0][0])
        doi = item.get("DOI", "")
        journal = item.get("container-title", [""])[0]
        bibtex = _item_to_bibtex(item)

        results.append({
            "title": title,
            "authors": authors,
            "year": year,
            "doi": doi,
            "journal": journal,
            "bibtex": bibtex,
        })
    return results


def validate_bibtex_file(bib_path: str) -> dict:
    """Validate a BibTeX file by checking DOIs against CrossRef.

    Args:
        bib_path: Path to a .bib file

    Returns:
        Dict with keys: total, valid, errors, warnings, corrected_entries
    """
    with open(bib_path, "r") as f:
        content = f.read()

    # Simple regex-based BibTeX parser (handles common cases)
    entries = re.findall(
        r"@(\w+)\{([^,]+),\s*(.*?)\n\}",
        content,
        re.DOTALL,
    )

    results = {
        "total": len(entries),
        "valid": 0,
        "errors": [],
        "warnings": [],
        "corrected_entries": [],
    }

    for entry_type, key, fields_str in entries:
        # Extract DOI if present
        doi_match = re.search(r"doi\s*=\s*\{([^}]+)\}", fields_str, re.IGNORECASE)
        title_match = re.search(r"title\s*=\s*\{([^}]+)\}", fields_str, re.IGNORECASE)

        if doi_match:
            doi = doi_match.group(1).strip()
            try:
                correct_bib = doi_to_bibtex(doi)
                results["valid"] += 1
                results["corrected_entries"].append({
                    "key": key.strip(),
                    "original_doi": doi,
                    "corrected_bibtex": correct_bib,
                })
            except ValueError as e:
                results["errors"].append({
                    "key": key.strip(),
                    "error": f"Invalid DOI: {e}",
                })
        elif title_match:
            title = title_match.group(1).strip()
            # Try to find the paper by title
            search_results = search_papers(title, rows=1)
            if search_results and "error" not in search_results[0]:
                best = search_results[0]
                results["warnings"].append({
                    "key": key.strip(),
                    "warning": "No DOI in entry — found possible match",
                    "suggested_doi": best["doi"],
                    "suggested_bibtex": best["bibtex"],
                })
                results["valid"] += 1
            else:
                results["errors"].append({
                    "key": key.strip(),
                    "error": "No DOI and could not find paper by title",
                })
        else:
            results["errors"].append({
                "key": key.strip(),
                "error": "No DOI or title found in entry",
            })

    return results


def find_and_cite(
    description: str,
    n: int = 3,
    from_year: int | None = None,
) -> str:
    """Search for papers matching a description and return BibTeX entries.

    This is the main entry point for the Trainee: describe what you need
    to cite and get back verified BibTeX entries.

    Args:
        description: What the citation should support (e.g., "CRISPR gene editing
                     efficiency in human cells")
        n: Number of citations to return
        from_year: Only return papers from this year onward

    Returns:
        A string containing n BibTeX entries, ready to paste into references.bib
    """
    results = search_papers(description, rows=n, filter_from_year=from_year)
    bibtex_entries = []
    for r in results:
        if "error" not in r and r.get("bibtex"):
            bibtex_entries.append(r["bibtex"])
    return "\n\n".join(bibtex_entries) if bibtex_entries else "% No papers found for this query"
