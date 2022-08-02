using NonNegLeastSquares
α = 1.3
β = 1.021389879988043
γ = 475.0377082498354
γ = sqrt(475.0377082498354)

A = [1. 1. 1. 1.; #iter0 
     α  β  γ  1.; #iter129
     α  1. γ  1.; #no forcing adjustment
     1. β  1. 1.#no adjustemnt to init or diff
     ]
ΔH = [1.0;
      round(2.526668450213555, digits = 2);
      round(2.117074817625146, digits = 2);
      round(2.418236154861419, digits = 2);
      ]

# coefs1 = A \ ΔH
coefs2 =nonneg_lsq(A, ΔH)

#using relative change to quantify these

α = 0.3
β = 0.04
γ = 0.4221480191380372
# γ = -0.09417824611240679 #airtemperature
A = [1. 1. 1. 1.; #iter0 
     α  β  γ  1.; #iter129
     α  1. γ  1.; #no forcing adjustment
     1. β  1. 1.#no adjustemnt to init or diff
     ]
ΔH = [1.0;
      round(0.8657850726631864, digits = 5);
      round(0.7167455919304684, digits = 5);
      round(0.8298058358808253, digits = 5);
      ]

coefs1 = A \ ΔH
coefs2 =nonneg_lsq(A, ΔH)
