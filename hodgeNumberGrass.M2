loadPackage "Schubert2"
G = flagBundle {2,2} -- secretly Gr(2,5)
X = sectionZeroLocus(OO_G(2)) -- quadric and a cubic
OmG = cotangentBundle G
OmX = cotangentBundle X

for i from 0 to floor(dim X / 2) list
   abs(chi exteriorPower(i, OmX) - chi exteriorPower(i, OmG))