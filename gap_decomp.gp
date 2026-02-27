# [ "A11", "A6", "A7", "L2(11)", "M11", "U4(2)", "U5(2)" ]
# G := SimpleGroup("A7"); yes
# G := SimpleGroup("L2(11)"); yes
G := SimpleGroup("A7"); 
# G := PSL(2,11);;



irrs := Irr(G);;
classes := ConjugacyClasses(G);;

Print("Irrep dimensions: ", List(irrs, chi -> chi[1]), "\n");;

V_irr := Filtered(irrs, chi -> chi[1] = 10)[2];;
Print("V character values: ", ValuesOfClassFunction(V_irr), "\n");;


ClassIndex := function(g)
    return First([1..Length(classes)], i -> g in classes[i]);
end;;

wedge3_vals := List([1..Length(classes)], function(i)
     local g;
     g := Representative(classes[i]);
     return (V_irr[i]^3
             - 3 * V_irr[ClassIndex(g^2)] * V_irr[i]
             + 2 * V_irr[ClassIndex(g^3)]) / 6;
 end);;

wedge3_char := ClassFunction(G, wedge3_vals);;

Print("Decomposition of wedge^3 V:\n");;

for irr in irrs do
     mult := ScalarProduct(wedge3_char, irr);;
     if mult <> 0 then
         Print("  dim ", irr[1], " with multiplicity ", mult, "\n");;
     fi;;
 od;;