
using Plots, BenchmarkTools, Polynomials
pyplot(size=(600,600))

#-------------------------------------------------------------------------------
function integral(f; h = 1e-2)
  g(x) = f(tanh(sinh(x)))*(cosh(x)*(sech(sinh(x))^2))
  i = 2*h
  S = g(0)
  while tanh(sinh(i)) != 1
    S += g(i) + g(-i)
    i += 2*h
  end
  return 2*h*S
end

function integral(f, a, b; h = 1e-2)
  c, δ = BigFloat(a + b) / 2, BigFloat(b - a) / 2
  g(t) = f(c + δ * t) * δ
  return integral(g, h=h)
end
#------------------------------------------------------------------------------
function coef_b(w, fi, i_1, i_2)
    cima = integral(x-> x*w(x)*(fi(x)^2), i_1, i_2)
    baixo = integral(x-> w(x)*((fi(x)^2)), i_1, i_2)
    b = cima/baixo
    return round(b;digits=8)
end

function coef_c(w, fi_1, fi_2, i_1, i_2)
    cima = integral(x-> x*w(x)*fi_1(x)*fi_2(x), i_1, i_2)
    baixo = integral(x-> w(x)*((fi_2(x)^2)), i_1, i_2)
    c = cima/baixo
    return round(c;digits=8)
end
#------------------------------------------------------------------------------
function Fpoliort(n, w, i_1, i_2)
    F = []
    B = []
    C = []
    p = Poly([1])
    push!(F,p)
    b = coef_b(w, F[1], i_1, i_2)
    push!(B, b)
    p = Poly([-b,1])
    push!(F, p)
    for k = 2:n
        b = coef_b(w, F[k], i_1, i_2)
        push!(B, b)
        c = coef_c(w, F[k], F[k-1], i_1, i_2)
        push!(C, c)
        p = (Poly([-B[k],1])*F[k]) - (C[k-1]*F[k-1])
        push!(F,p)
    end
    return F
end
#-------------------------------------------------------------------------------
function Apoliort(F,f,n,i_1,i_2)
    M = zeros(n+1,n+1)
    A = zeros(n+1)
    B = zeros(n+1)
    for j = 0:n
            M[j+1,j+1] = integral(x -> w(x)*F[j+1](x)*F[j+1](x),i_1,i_2)
    end
    for j = 0:n
        B[j+1] = integral(x -> w(x)*f(x)*F[j+1](x),i_1,i_2)
    end
    A = M\B
    return round.(A;digits=6)
end
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
function Ppoliort(A,F)
    p = Poly([0])
    m = length(A)
    for i = 1:m
        p = p + A[i]*F[i]
    end
    coef = round.(Float64.(coeffs(p));digits=6)
    p = Poly(coef)
    return p
end
#------------------------------------------------------------------------------
function Apoliquad(f,n,i_1,i_2)
    M = zeros(n+1,n+1)
    A = zeros(n+1)
    B = zeros(n+1)
    for j = 0:n
        for k = j:n
            M[j+1,k+1] = integral(x -> x^(j+k),i_1,i_2)
        end
    end
    for j = 1:n
        for k = 0:j-1
            M[j+1,k+1] = M[k+1,j+1]
        end
    end
    for j = 0:n
        B[j+1] = integral(x -> (x^(j))*f(x),i_1,i_2)
    end
    A = M\B
    return round.(A;digits=6)
end
#------------------------------------------------------------------------------
function Vpoliquad(x,A)
    v = 0.0
    m = length(A)
    for i = 1:m
        v += A[i]*x^(i-1)
    end
    return v
end
#------------------------------------------------------------------------------
#Entradas:
n = 30
f = x-> log(x)
w = x-> 1/sqrt(1-x^2)
i_1 = 0
i_2 = 1

F_1=Fpoliort(n, w, i_1, i_2)
A_1=Apoliort(F_1,f,n,i_1,i_2)
A_2=Apoliquad(f,n,i_1,i_2)

p = Ppoliort(A_1,F_1)


#Execução principal
F_1 = Fpoliort(n, w, i_1, i_2)
A_1 = Apoliort(F_1,f,n,i_1,i_2)
A_2 = Apoliquad(f,n,i_1,i_2)

#Plot da função e dos polinômios
plot(x->f(x), i_1, i_2, c=:black, label=:"fun")
println(@benchmark plot!(x->p(x), i_1, i_2, c=:red, label=:"ort"))
println(@benchmark plot!(x->Vpoliquad(x,A_2), i_1, i_2, c=:blue, label=:"quad"))






