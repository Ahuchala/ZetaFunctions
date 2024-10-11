f = "x01^2+2*x02^2+4*x03^2+5*x04^2+6*x12^2+11*x13^2+75*x14^2+13*x23^2+43*x24^2+8*x34^2"

n = 5


import re

def replace_index_in_polynomial(polynomial, i, j):
    # Define the regular expression to match x followed by multi-index
    pattern = re.compile(r"(x(\d+))")
    
    def replace_index(match):
        # Extract the multi-index part (group 2)
        multi_index = match.group(2)
        
        # Replace i with j in the multi-index
        new_index = multi_index.replace(str(i), str(j))
        
        # Check if there are repeated digits in the new index
        if len(new_index) != len(set(new_index)):
            return "0"  # Replace the entire monomial with 0 if repeated digits exist
        
        # Return the updated variable
        return f"x{new_index}"
    
    # Apply the replacement on the entire polynomial
    updated_polynomial = pattern.sub(replace_index, polynomial)
    
    # Clean up any monomials that became zero
    updated_polynomial = re.sub(r"\b0[\^0-9]*\b", "0", updated_polynomial)
    
    return updated_polynomial

# Example usage
i = 3  
j = 4


ls = []
for i in range(n):
    for j in range(n):
        if i != j:
            ls.append(replace_index_in_polynomial(f,i,j))
            ls.append(replace_index_in_polynomial(f,i,i) + "-"+replace_index_in_polynomial(f,j,j))

print(ls)
