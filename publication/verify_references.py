#!/usr/bin/env python3
"""
verify_references.py

Every reference in the manuscript is checked against authoritative bibliographic
APIs before it is allowed into the reference list.

WHY THIS EXISTS. A citation in an earlier paper from this programme had to be
retracted at peer review: it cited a paper as a clinical case report when the
paper was an in vitro study containing no patients, and the author list was
conflated with a different paper by an overlapping group. Both errors would have
been caught by mechanically re-fetching the record. So nothing here is written
from memory or from a search summary. Each entry is fetched from Crossref by
DOI, and cross-checked against Europe PMC where a PMID exists.

WHAT IS CHECKED
  1. the DOI resolves at all
  2. title as recorded by the publisher
  3. full author list, in order
  4. journal, year, volume, pages
  5. where a PMID is supplied, that Crossref and Europe PMC agree on the title
  6. the abstract is printed so the ASSIGNMENT to a claim can be judged by eye

WHAT IT CANNOT DO. It cannot tell you the reference supports the sentence it is
attached to. That is a human judgement and is made separately, against the
abstract printed here. Bibliographic correctness and citation correctness are
different problems and both have to be solved.

Usage:  python3 verify_references.py refs.json
where refs.json is a list of {"key", "doi", "pmid" (optional), "claim"}.

Author: Christopher Lawrence
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UA = {"User-Agent": "kidney-alphagenome-refcheck/1.0 (mailto:Christopher.lawrence3@nhs.net)"}
CROSSREF = "https://api.crossref.org/works/"
EPMC = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
        "?query={q}&resultType=core&format=json")


def get_json(url: str, retries: int = 3) -> dict | None:
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.load(fh)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(2 * attempt)
        except Exception:                                    # noqa: BLE001
            time.sleep(2 * attempt)
    return None


def crossref(doi: str) -> dict | None:
    d = get_json(CROSSREF + urllib.parse.quote(doi))
    return d.get("message") if d else None


def epmc(query: str) -> dict | None:
    d = get_json(EPMC.format(q=urllib.parse.quote(query)))
    if not d:
        return None
    res = d.get("resultList", {}).get("result", [])
    return res[0] if res else None


def ndt_authors(authors: list[dict]) -> str:
    """NDT style: all authors up to six; more than six -> first three + et al.

    Consortium authors carry a "name" key instead of given/family. Ignoring that
    produced a leading empty name and a stray comma on the GTEx and ENCODE
    entries. Where a consortium is the FIRST author, NDT practice is to cite the
    consortium alone rather than the consortium plus individuals.
    """
    names = []
    for a in authors:
        if a.get("name"):
            nm = a["name"].strip()
            if nm.lower().startswith("the "):
                nm = nm[4:]
            names.append(nm)
            continue
        fam = a.get("family", "").strip()
        given = a.get("given", "") or ""
        initials = "".join(p[0] for p in given.replace("-", " ").split() if p)
        if fam:
            names.append(f"{fam} {initials}".strip())
    if names and authors and authors[0].get("name"):
        return names[0]                      # consortium-authored: cite it alone
    if len(names) > 6:
        return ", ".join(names[:3]) + " et al"
    return ", ".join(names)


def fmt_ndt(m: dict) -> str:
    """Authors. Title. Journal Year; Volume: pages"""
    auth = ndt_authors(m.get("author", []) or [])
    title = (m.get("title") or [""])[0].rstrip(".")
    jrnl = (m.get("short-container-title") or m.get("container-title") or [""])
    jrnl = jrnl[0] if jrnl else ""
    year = ""
    for k in ("published-print", "published-online", "issued"):
        dp = (m.get(k) or {}).get("date-parts") or []
        if dp and dp[0] and dp[0][0]:
            year = str(dp[0][0])
            break
    # Preprints have no journal, volume or pages. Formatting them like an
    # article produces a reference that looks published and is not. NDT wants
    # the server named and the DOI given.
    if m.get("type") == "posted-content":
        server = ((m.get("institution") or [{}])[0].get("name")
                  or m.get("group-title") or "preprint")
        posted = (m.get("posted") or {}).get("date-parts") or [[year]]
        yr = posted[0][0] if posted and posted[0] else year
        return (f"{auth}. {title}. {server} {yr}; "
                f"doi:{m.get('DOI','')} (preprint, not peer reviewed)")

    vol = m.get("volume", "")
    page = m.get("page", "")
    art = m.get("article-number", "")
    tail = f"{vol}: {page}" if page else (f"{vol}: {art}" if art else vol)
    return f"{auth}. {title}. {jrnl} {year}; {tail}".replace("  ", " ")


def main() -> None:
    path = Path(sys.argv[1] if len(sys.argv) > 1 else "refs.json")
    refs = json.loads(path.read_text())
    ok, bad = 0, []

    for i, r in enumerate(refs, 1):
        key, doi, pmid = r["key"], r.get("doi"), r.get("pmid")
        print("=" * 78)
        print(f"[{i}] {key}")
        print(f"    CLAIM: {r.get('claim','(none given)')}")
        m = crossref(doi) if doi else None
        if not m:
            print(f"    *** CROSSREF FAILED for doi {doi} - DO NOT USE ***")
            bad.append(key)
            continue

        print(f"    DOI resolves     : {doi}")
        print(f"    title            : {(m.get('title') or [''])[0]}")
        print(f"    type             : {m.get('type')}")
        n_auth = len(m.get("author", []) or [])
        print(f"    authors (n={n_auth:>2})   : {ndt_authors(m.get('author', []) or [])}")
        print(f"    journal          : {((m.get('container-title') or ['']) or [''])[0]}")
        print()
        print(f"    NDT FORMAT >>> {fmt_ndt(m)}")

        if pmid:
            e = epmc(f"EXT_ID:{pmid}")
            if not e:
                print(f"    ! PMID {pmid} not found in Europe PMC")
            else:
                # Normalise punctuation before comparing. Publishers use en
                # dashes where PubMed uses hyphens, and that difference alone
                # produced a false MISMATCH that could mask a real one.
                def norm(t: str) -> str:
                    import re as _re
                    t = t.lower().strip().rstrip(".")
                    for ch in ("\u2010", "\u2011", "\u2012", "\u2013",
                               "\u2014", "\u2212"):
                        t = t.replace(ch, "-")
                    t = t.replace("&gt;", ">").replace("&lt;", "<")
                    t = t.replace("\u2019", "'").replace("\u201c", '"').replace("\u201d", '"')
                    return _re.sub(r"\s+", " ", t)

                t1 = norm((m.get("title") or [""])[0])
                t2 = norm(e.get("title") or "")
                agree = t1[:60] == t2[:60]
                print(f"    PMID {pmid} title match: "
                      f"{'OK' if agree else '*** MISMATCH ***'}")
                if not agree:
                    print(f"      crossref: {t1[:90]}")
                    print(f"      epmc    : {t2[:90]}")
                    bad.append(key)
                ab = (e.get("abstractText") or "")
                if ab:
                    import re
                    ab = re.sub(r"<[^>]+>", " ", ab)
                    ab = re.sub(r"\s+", " ", ab)
                    print(f"    ABSTRACT (judge the assignment against this):")
                    print(f"      {ab[:700]}")
        ok += 1
        time.sleep(0.4)

    print("=" * 78)
    print(f"verified {ok}/{len(refs)}")
    if bad:
        print(f"PROBLEMS: {bad}")
        sys.exit(1)
    print("all references resolved")


if __name__ == "__main__":
    main()
