# Always-over-C irreps via Dixon, filtered by dimension, printing generator images.

PrintComplexIrrepsByDim := function(G, dim)
  local gens, reps, good, k, rho, i, M, row, deg;

  if not IsGroup(G) then
    Error("Input must be a GAP group object, not a character table name.");
  fi;

  gens := GeneratorsOfGroup(G);
  if Length(gens) = 0 then
    Error("Group has no generators?");
  fi;

  # Dixon irreducibles over cyclotomics (conceptually "over C")
  reps := IrreducibleRepresentationsDixon(G);

  deg := r -> Length(Image(r, gens[1]));  # matrix size = representation degree
  good := Filtered(reps, r -> deg(r) = dim);

  Print("Found ", Length(good), " irreps of degree ", dim, ".\n\n");

  for k in [1..Length(good)] do
    rho := good[k];
    Print("=== irrep ", k, " (degree ", dim, ") ===\n");
    for i in [1..Length(gens)] do
      Print("Matrix image of generator ", i, ":\n");
      M := Image(rho, gens[i]);
      for row in M do
        Print(row, "\n");
      od;
      Print("\n");
    od;
  od;

  return good;
end;

G := AlternatingGroup(7);;
Size(G);   # 720

PrintComplexIrrepsByDim(G, 10);   # prints the 10-dimensional complex irreps (if any)