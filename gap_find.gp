#############################################################################
# Search character tables for degree-n complex representations V
# such that wedge^k(V) has exactly ONE 1-dimensional constituent (any linear).
# Prints V and wedge^k(V) decompositions.
#
# Default use-case: n=10, k=3
#############################################################################

LoadPackage("ctbllib");

#############################################################################
# PARAMETERS (edit these)
#############################################################################
PARAM_n := 10;              # dimension of V
PARAM_k := 3;               # exterior power
MAX_GROUP_SIZE := 5000;     # only test tables with |G| < this
MAX_CLASSES := 200;         # only test tables with <= this many classes
PRINT_WEDGE_DECOMP := true; # print full wedge^k decomposition
PRINT_PROGRESS_EVERY := 20; # set 0 to disable

#############################################################################
# Helpers
#############################################################################

MainId := function(tbl)
  if IsDuplicateTable(tbl) then
    return IdentifierOfMainTable(tbl);
  else
    return Identifier(tbl);
  fi;
end;

DecompByScalarProducts := function(tbl, phi)
  local irrs, i, m, out;
  irrs := Irr(tbl);
  out := [];
  for i in [1..Length(irrs)] do
    m := ScalarProduct(tbl, phi, irrs[i]);
    if m <> 0 then
      Add(out, rec(i := i, deg := irrs[i][1], mult := m));
    fi;
  od;
  return out;
end;

PrintDecomp := function(decomp)
  local j;
  if Length(decomp) = 0 then
    Print("0");
    return;
  fi;
  for j in [1..Length(decomp)] do
    if j > 1 then Print(" + "); fi;
    Print(decomp[j].mult, "*chi_", decomp[j].i, "[deg ", decomp[j].deg, "]");
  od;
end;

Total1DimMultiplicity := function(tbl, phi)
  local irrs, i, tot;
  irrs := Irr(tbl);
  tot := 0;
  for i in [1..Length(irrs)] do
    if irrs[i][1] = 1 then
      tot := tot + ScalarProduct(tbl, phi, irrs[i]);
    fi;
  od;
  return tot;
end;

# Build a virtual character chi from a witness record (keep/mults)
BuildCharacterFromWitness := function(tbl, wit)
  local irrs, chi, j;
  irrs := Irr(tbl);
  chi := 0 * irrs[1];
  for j in [1..Length(wit.keep)] do
    if wit.mults[j] <> 0 then
      chi := chi + wit.mults[j] * irrs[ wit.keep[j] ];
    fi;
  od;
  return chi;
end;

#############################################################################
# Core: enumerate ALL degree-n reps (as nonnegative sums of irreds)
# using only irreps with degree <= n, and test wedge^k condition.
#
# This is streaming/backtracking: no huge list of solutions stored unless found.
#############################################################################

FindAllWitnessesDegreeN_WedgeK_One1D := function(tbl, n, k)
  local irrs, degs, keep, irrs2, degs2, mmax, mults, chi, ans, recfun, d, m, oldchi, wk;

  irrs := Irr(tbl);
  degs := List(irrs, x -> x[1]);

  # irreps of degree > n cannot appear in a degree-n sum
  keep := Filtered([1..Length(irrs)], function(i) return degs[i] <= n; end);
  if Length(keep) = 0 then
    return [];
  fi;

  irrs2 := List(keep, i -> irrs[i]);
  degs2 := List(keep, i -> degs[i]);

  mults := List([1..Length(irrs2)], i -> 0);
  chi := 0 * irrs2[1];
  ans := [];

  recfun := function(pos, remaining)
    if remaining = 0 then
      # chi(1) should be n automatically; keep a check anyway
      if chi[1] = n then
        wk := ExteriorPower(chi, k);
        if Total1DimMultiplicity(tbl, wk) = 1 then
          Add(ans, rec(keep := ShallowCopy(keep), mults := ShallowCopy(mults)));
        fi;
      fi;
      return;
    fi;

    if pos > Length(irrs2) then
      return;
    fi;

    d := degs2[pos];
    mmax := Int(remaining / d);

    for m in [0..mmax] do
      mults[pos] := m;

      oldchi := chi;
      if m <> 0 then
        chi := chi + m * irrs2[pos];
      fi;

      recfun(pos+1, remaining - m*d);

      chi := oldchi;
      mults[pos] := 0;
    od;
  end;

  recfun(1, n);
  return ans;
end;

#############################################################################
# Main scan: collect unique tables within bounds, sort small-first, scan.
#############################################################################

ScanAndPrintSolutions := function(n, k, maxSize, maxClasses)
  local names, seen, uniq, id, tbl, key, i, t, wits, s, chi, wk;

  names := AllCharacterTableNames();
  seen := [];
  uniq := [];

  # collect candidate tables (dedupe + bounds)
  for id in names do
    tbl := CharacterTable(id);

    if Size(tbl) < maxSize and NrConjugacyClasses(tbl) <= maxClasses then
      key := MainId(tbl);
      if not key in seen then
        Add(seen, key);
        Add(uniq, tbl);
      fi;
    fi;
  od;

  # sort: |G| then #classes then id
  Sort(uniq, function(a, b)
    if Size(a) <> Size(b) then
      return Size(a) < Size(b);
    elif NrConjugacyClasses(a) <> NrConjugacyClasses(b) then
      return NrConjugacyClasses(a) < NrConjugacyClasses(b);
    else
      return MainId(a) < MainId(b);
    fi;
  end);

  Print("Testing ", Length(uniq), " unique tables (sorted by |G|)...\n");
  Print("Target: dim(V) = ", n, ", wedge^", k, " V has exactly one 1-d constituent.\n\n");

  for i in [1..Length(uniq)] do
    t := uniq[i];

    if PRINT_PROGRESS_EVERY > 0 and i mod PRINT_PROGRESS_EVERY = 0 then
      Print("... progress ", i, "/", Length(uniq),
            " (currently ", MainId(t),
            ", |G|=", Size(t),
            ", classes=", NrConjugacyClasses(t), ")\n");
    fi;

    wits := FindAllWitnessesDegreeN_WedgeK_One1D(t, n, k);

    if Length(wits) > 0 then
      Print("HIT: ", MainId(t),
            " |G|=", Size(t),
            " classes=", NrConjugacyClasses(t),
            "  (#solutions=", Length(wits), ")\n");

      for s in [1..Length(wits)] do
        chi := BuildCharacterFromWitness(t, wits[s]);
        wk  := ExteriorPower(chi, k);

        Print("  sol ", s, ": V = ");
        PrintDecomp(DecompByScalarProducts(t, chi));
        Print("   (total dim = ", chi[1], ")\n");

        Print("        total 1-d in wedge^", k, " V = ",
              Total1DimMultiplicity(t, wk), "\n");

        if PRINT_WEDGE_DECOMP then
          Print("        wedge^", k, " V = ");
          PrintDecomp(DecompByScalarProducts(t, wk));
          Print("\n");
        fi;
      od;

      Print("\n");
    fi;
  od;
end;

#############################################################################
# RUN
#############################################################################

ScanAndPrintSolutions(PARAM_n, PARAM_k, MAX_GROUP_SIZE, MAX_CLASSES);