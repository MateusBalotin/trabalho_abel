# trabalho_abel
using Plots
pyplot(size=(800,800))

#Método de integração (posteriormente implementar tansin, comparar qual é melhor)
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


#-------------------------------------------------------------------------------
#Funções coeficientes para polinômio ortogonal
function coef_a(w, fi, f, i_1, i_2)
    cima = integral(x-> w(x)*fi(x)*f(x), i_1, i_2)
    baixo = integral(x-> w(x)*((fi(x)^2)), i_1, i_2)
    a = cima/baixo
    return a
end

function coef_b(w, fi, i_1, i_2)
    cima = integral(x-> x*w(x)*(fi(x)^2), i_1, i_2)
    baixo = integral(x-> w(x)*((fi(x)^2)), i_1, i_2)
    b = cima/baixo
    return b
end

function coef_c(w, fi_1, fi_2, i_1, i_2)
    cima = integral(x-> x*w(x)*fi_1(x)*fi_2(x), i_1, i_2)
    baixo = integral(x-> w(x)*((fi_2(x)^2)), i_1, i_2)
    c = cima/baixo
    return c
end
#Definindo o polinômio
function poliort(n, w, f, i_1, i_2)
    A = []
    B = []
    C = []
    F = Any[x-> 1]
    #Calculando ao, b1, fi1, a1
    a = coef_a(w, F[1], f, i_1, i_2)
    push!(A, a)
    b = coef_b(w, F[1], i_1, i_2)
    push!(B, b)
    fi = x -> x - b
    push!(F, fi)
    a = coef_a(w, F[2], f, i_1, i_2)
    push!(A, a)
    #Calculando ai's e fii's
    for k = 2:n
        b = coef_b(w, F[k], i_1, i_2)
        push!(B, b)
        c = coef_c(w, F[k], F[k-1], i_1, i_2)
        push!(C, c)
        fih = x -> (x-B[k])*F[k](x) - (C[k-1]*F[k-1](x))
        push!(F,fih)
        a = coef_a(w, F[k+1], f, i_1, i_2)
        push!(A, a)
    end
    return A, F
end
#Aplicando em x
function vpoliort(x)
    v = 0.0
    A,F = poliort(n, w, f, i_1, i_2)
    m = length(A)
    for i = 1:m
        v += A[i]*F[i](x)
    end
    return v
end
#-------------------------------------------------------------------------------
#Quadrados Mínimos
#Definindo o polinômio
function poliquad(f,n,i_1,i_2)
    I = ones(n,n)
    A = ones(n)
    B = ones(n)
    for j = 1:n
        for k = 1:n
            I[j,k] = integral(x -> x^(j+k),i_1,i_2)
        end
    end
    for j = 1:n
        B[j] = integral(x -> (x^j)*f(x),i_1,i_2)
    end
    A = I\B
    return A
end
#Aplicando em x
function vpoliquad(x)
    v = 0.0
    A = poliquad(f,n,i_1,i_2)
    for i = 1:n
        v += A[i]*x^(i)
    end
    return v
end
#-------------------------------------------------------------------------------
function erro()
    erro_ort = integral(x->(f(x)-vpoliort(x))^2, i_1, i_2)
    erro_quad = integral(x->(f(x)-vpoliquad(x))^2, i_1, i_2)
end
#------------------------------------------------------------------------------    
n = 6
i_1 = -1
i_2 = 1
f = x -> 1/sqrt(1-x^2)
w = x -> 1/sqrt(1-x^2)

plot(x->vpoliquad(x),i_1,i_2,c=:blue)
plot!(x->vpoliort(x),i_1,i_2,c=:black) 
plot!(x->f(x),i_1,i_2,c=:green) 
