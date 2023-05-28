function f0(m)
    return sum([Int(floor(log(max(1,m-i))/log(p))) for i = 1 :n])
end

# be careful you're passing in a bigint
function zp_val(m)
    ans = 0
    while mod(m,p) == 0
        ans += 1
        m ÷= p
    end
    return ans
end

function g_prec(m,i)
    return zp_val(binomial(big(-m),big(i))) 
end

# this might return too long of a list...
# integer m, list A0
function make_A(m,A0)
    j,j1,N,l = -1,-1,-1,-1
    A = Dict()
    Aprime = Dict()
    for i = 1:max(m,n)
        if haskey(A0,i)
            A[i] = A0[i]
        else
            A[i] = f0(i)
        end
    end
    go_to_step_2 = true
    go_to_step_3 = false
    go_to_step_4 = false
    go_to_step_5 = false
    go_to_step_6 = false
    go_to_step_7 = false

    while true
        if go_to_step_2
            Aprime = A
            j = n+1
            go_to_step_2 = false
            go_to_step_3 = true
        end
        if go_to_step_3
            if !haskey(A,j)
                if A == Aprime
                    return A
                end
            end
            go_to_step_3 = false
            go_to_step_4 = true
        end
        if go_to_step_4
            j1 = p*Int(ceil(j/p))
            N = n - 1 + A[j1÷p]
            l = 1
            go_to_step_4 = false
            go_to_step_5 = true
        end
        if go_to_step_5
            if n*p< j1 +l*p && n*log(j1+l*n)/log(p)-l <= N
                for blah = j:j1
                    A[blah] = min(N,f0(j1))
                end
                j = j1+1
                go_to_step_5 = false
                go_to_step_3 = true
            else
                go_to_step_5 = false
                go_to_step_6 = true
            end
        end
        if go_to_step_6
            if f0(j1+l*p)-l-g_prec(j1,l)<=N
                l = l+1
                go_to_step_6 = false
                go_to_step_5 = true
            else
                go_to_step_6 = false
                go_to_step_7 = true
            end
        end
        if go_to_step_7
            len_A_keys = size(collect(keys(A)),1)
            for i = len_A_keys+1:j1+l*p
                A[i] = f0(i)
            end
            N = max(A[j1+l*p]-l-g_prec(j1,l),N)
            l = l+1
            go_to_step_7 = false
            go_to_step_5 = true
        end
    end
end


function find_s(r)
    s = r
    A = Dict()
    j = -1
    go_to_step_2 = true
    go_to_step_3 = false

    while true
        if go_to_step_2
            j = s - n + 1
            go_to_step_2 = false
            go_to_step_3 = true
        end
        if go_to_step_3
            if j > 0 && n*log(p*(n+j)-1)/log(p)<=n-1+j-r
                return s
            end
            go_to_step_3 = false
            go_to_step_4 = true
        end
        if go_to_step_4
            A = make_A(p*(n+j),A)
            if A[p*(n+j)]>n-1+j-r
                s = s+1
                go_to_step_4 = false
                go_to_step_2 = true
            else
                j = j+1
                go_to_step_4 = false
                go_to_step_3 = true
            end
        end
    end
end


