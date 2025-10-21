#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate a reference Newick tree with bootstrap supports computed
from a file of bootstrap replicate trees (e.g., RAxML[-NG] *.bootstraps).
- Keeps the topology and branch lengths of the reference tree intact.
- Replaces internal node labels with integer bootstrap percentages.
- Robust to minor naming issues in bootstrap trees:
    * Duplicate leaf names (e.g., "G-IBV1" appearing twice) are
      automatically mapped to likely missing counterparts (e.g., "G-IBV2").
    * A few stray or mismatched labels can be aliased via an optional TSV file.
Usage
-----
python tree_apply_bootstrap.py \
  --ref treerooted.nwk \
  --boot RAxML_tree.raxml.bootstraps \
  --out treerooted_with_bootstrap.nwk \
  --log bootstrap_mapping.log
Optional:
  --aliases aliases.tsv   # two columns: <bootstrap_name>\t<reference_name>
Notes
-----
If you prefer the official mapping, you can also use:
  raxml-ng --support --tree treerooted.nwk --bs-trees RAxML_tree.raxml.bootstraps
This script exists to fix difficult label issues and to keep everything in-Python.
"""
from __future__ import annotations
import argparse, re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional, List, Set, Dict

# --------- Minimal Newick utilities (no external deps) ---------
class Node:
    __slots__ = ("name","length","children","parent")
    def __init__(self, name: Optional[str]=None, length: Optional[str]=None):
        self.name = name
        self.length = length
        self.children: List['Node'] = []
        self.parent: Optional['Node'] = None
    def add_child(self, child: 'Node'):
        self.children.append(child); child.parent = self
    def is_leaf(self)->bool:
        return len(self.children)==0

def _parse_newick(s:str)->Node:
    s = re.sub(r"\s+","", s.strip())
    if not s.endswith(";"):
        raise ValueError("Newick string must end with ';'")
    n = len(s); i = 0; stack: List[Node] = []; root: Optional[Node] = None
    def read_name(j:int):
        start = j
        while j < n and s[j] not in ",:()":
            if s[j] == ';': break
            j += 1
        return s[start:j], j
    def read_len(j:int):
        if j < n and s[j] == ':':
            j += 1; start = j
            while j < n and s[j] not in ",()":
                if s[j] == ';': break
                j += 1
            return s[start:j], j
        return None, j
    while i < n:
        c = s[i]
        if c == '(':
            stack.append(Node()); i += 1
        elif c == ',':
            i += 1
        elif c == ')':
            node = stack.pop(); i += 1
            name, i = read_name(i)
            if name: node.name = name
            length, i = read_len(i)
            if length: node.length = length
            if stack: stack[-1].add_child(node)
            else: root = node
        elif c == ';':
            i += 1; break
        else:
            name, i = read_name(i)
            length, i = read_len(i)
            leaf = Node(name=name or None, length=length)
            if stack: stack[-1].add_child(leaf)
            else: root = leaf
    if root is None: raise ValueError("Failed to parse Newick")
    return root

def parse_newick(nwk_text:str)->Node:
    return _parse_newick(nwk_text)

def _leaf_list(n:Node)->List[str]:
    if n.is_leaf(): return [n.name]
    out = []
    for ch in n.children: out += _leaf_list(ch)
    return out

def get_leaf_labels(n:Node)->Set[str]:
    return set(_leaf_list(n))

def _tipset(n:Node)->Set[str]:
    if n.is_leaf(): return {n.name}
    s = set()
    for ch in n.children: s |= _tipset(ch)
    return s

def canonical_split(a:Set[str], all_tips:Set[str]):
    b = all_tips - a
    if 0 < len(a) < len(all_tips):
        # canonical: smaller by size (then lexicographically)
        if len(a) < len(b) or (len(a)==len(b) and tuple(sorted(a))<=tuple(sorted(b))):
            return frozenset(a)
        else:
            return frozenset(b)
    return frozenset()

def edge_splits(root:Node, all_tips:Set[str]):
    splits = set()
    def post(n:Node)->Set[str]:
        if n.is_leaf(): return {n.name}
        tips = set()
        for ch in n.children:
            t = post(ch); tips |= t
            cs = canonical_split(t, all_tips)
            if cs: splits.add(cs)
        return tips
    post(root)
    return splits

def write_newick(n:Node)->str:
    def rec(x:Node)->str:
        if x.is_leaf():
            label = x.name or ""
            return f"{label}:{x.length}" if x.length is not None else label
        inner = ",".join(rec(c) for c in x.children)
        label = x.name or ""
        if x.length is not None:
            return f"({inner}){label}:{x.length}"
        else:
            return f"({inner}){label}"
    return rec(n) + ";"

# --------- Bootstrap support computation ---------
def split_base_num(name:str):
    m = re.match(r'^(.*?)(\d+)$', name or "")
    return (m.group(1), int(m.group(2))) if m else (name, None)

def build_alias_map(ref_tips:Set[str], leaves:List[str], user_alias: Optional[Dict[str,str]]=None):
    counts = Counter(leaves)
    dups = [nm for nm,c in counts.items() if c>1]
    uniq = set(leaves)
    missing = sorted(ref_tips - uniq)
    extra = sorted(uniq - ref_tips)

    alias: Dict[str,str] = {}
    miss_set = set(missing)

    # Use user-provided aliases first
    if user_alias:
        alias.update(user_alias)
        # consider those resolved
        miss_set -= set(user_alias.values())

    # Resolve duplicates to likely missing counterparts (e.g., ...1 -> ...2)
    for d in dups:
        base, num = split_base_num(d)
        if num is not None:
            cands = [m for m in list(miss_set) if split_base_num(m)[0]==base and m!=d]
            if cands:
                tgt = sorted(cands)[0]
                alias[d+"#dup"] = tgt; miss_set.discard(tgt); continue
            for off in (1,-1,2,-2):
                cand = f"{base}{num+off}"
                if cand in miss_set:
                    alias[d+'#dup'] = cand; miss_set.discard(cand); break
        else:
            pref = (d.split('-')[0] if '-' in d else d[:1])
            cands = [m for m in list(miss_set) if (m.split('-')[0] if '-' in m else m[:1])==pref]
            if cands:
                tgt = sorted(cands)[0]
                alias[d+"#dup"] = tgt; miss_set.discard(tgt)

    # Map extras to remaining missing (try to keep same prefix letter)
    for ex in extra:
        if ex in alias: continue
        pref = (ex.split('-')[0] if '-' in ex else ex[:1])
        cands = [m for m in list(miss_set) if (m.split('-')[0] if '-' in m else m[:1])==pref]
        tgt = (sorted(cands)[0] if cands else (sorted(miss_set)[0] if miss_set else None))
        if tgt:
            alias[ex] = tgt; miss_set.discard(tgt)

    return alias, missing, extra, dups

def relabel_tree_leaves(root:Node, alias_map:Dict[str,str], leaf_counts:Counter):
    seen = defaultdict(int)
    def rec(n:Node):
        if n.is_leaf():
            nm = n.name; seen[nm]+=1
            key = nm
            if leaf_counts[nm]>1 and seen[nm]>1:
                key = nm+"#dup"
            n.name = alias_map.get(key, nm)
        else:
            for ch in n.children: rec(ch)
    rec(root)

def compute_support_map(boot_newicks:List[str], ref_tips:Set[str], log: List[str], user_alias: Optional[Dict[str,str]]=None):
    counts = Counter()
    used = 0
    for i, nwk in enumerate(boot_newicks):
        try:
            r = parse_newick(nwk)
        except Exception as e:
            log.append(f"[{i}] parse error: {e}"); continue
        leaves = _leaf_list(r); lc = Counter(leaves)
        uniq = set(leaves)
        if uniq != ref_tips or any(c>1 for c in lc.values()):
            alias, missing, extra, dups = build_alias_map(ref_tips, leaves, user_alias)
            relabel_tree_leaves(r, alias, lc)
            uniq2 = set(_leaf_list(r))
            if uniq2 != ref_tips:
                log.append(f"[{i}] mismatch after alias; missing={sorted(ref_tips-uniq2)}, extra={sorted(uniq2-ref_tips)}; alias={alias}")
                continue
            else:
                log.append(f"[{i}] resolved dups={dups}, extra={extra} with alias={alias}")
        splits = edge_splits(r, ref_tips)
        for sp in splits: counts[sp]+=1
        used += 1
    support = {sp: counts[sp]/used for sp in counts}
    return support, used

def annotate_support_on_reference(ref_root:Node, support_map:Dict, ref_tips:Set[str]):
    def rec(n:Node):
        for ch in n.children:
            tips_ch = _tipset(ch)
            sp = canonical_split(tips_ch, ref_tips)
            supp = support_map.get(sp, 0.0)
            if not ch.is_leaf():
                ch.name = str(int(round(100*supp)))
            rec(ch)
    rec(ref_root)

# --------- CLI ---------
def read_bootstrap_trees(path:Path)->List[str]:
    txt = path.read_text()
    trees = [t.strip()+";" for t in txt.split(";") if t.strip()]
    return trees

def read_aliases(path: Optional[Path]):
    if not path: return None
    alias = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"): continue
        parts = line.strip().split("\t")
        if len(parts)>=2:
            alias[parts[0]] = parts[1]
    return alias or None

def main():
    ap = argparse.ArgumentParser(description="Annotate a Newick tree with bootstrap supports.")
    ap.add_argument("--ref", required=True, help="Reference tree (Newick) whose topology/lengths will be kept.")
    ap.add_argument("--boot", required=True, help="Bootstrap trees file (multiple Newicks separated by ';').")
    ap.add_argument("--out", required=True, help="Output Newick with support values as internal node labels.")
    ap.add_argument("--log", default="bootstrap_mapping.log", help="Write mapping details here.")
    ap.add_argument("--aliases", default=None, help="Optional TSV with manual aliases: bootstrap_name<TAB>reference_name")
    args = ap.parse_args()

    ref_text = Path(args.ref).read_text()
    ref_root = parse_newick(ref_text)
    ref_tips = get_leaf_labels(ref_root)

    user_alias = read_aliases(Path(args.aliases)) if args.aliases else None

    boot_trees = read_bootstrap_trees(Path(args.boot))
    log: List[str] = []
    support_map, used = compute_support_map(boot_trees, ref_tips, log, user_alias)
    annotate_support_on_reference(ref_root, support_map, ref_tips)
    out_newick = write_newick(ref_root)
    Path(args.out).write_text(out_newick)
    Path(args.log).write_text("\n".join(log + [f"---", f"Used bootstrap trees: {used}/{len(boot_trees)}"]))
    print(f"Wrote: {args.out} (support from {used} bootstrap trees)")
    print(f"Log:   {args.log}")

if __name__ == "__main__":
    main()
