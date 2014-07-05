;lap_coeff.pro
;	by Joe Hahn, jhahn@spacescience.org

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Elliptic_K, x
;
;  Compute the complete elliptic integral of the first kind K(x), valid for
;  0 <= x < 1. From Abromowitz & Stegun. This is the integral of the function
;  f(x) = 1/sqrt( 1 - (x*sin(t))^2 ) from t=0 to t=Pi/2. Fractional errors are
;  smaller than 10^(-8) when compared to Maple.
;
m = x^2
m1 = 1d - m
K = 1.38629436112d        + 0.09666344259d*m1     + 0.03590092383d*(m1^2) + $
    0.03742563713d*(m1^3) + 0.01451196212d*(m1^4) +                         $
  ( 0.5d                  + 0.12498593597d*m1     + 0.06880248576d*(m1^2) + $
    0.03328355346d*(m1^3) + 0.00441787012d*(m1^4) )*alog(1d/m1)

return, K
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function Elliptic_E, x
;
;  Compute the complete elliptic integral of the second kind E(x), valid for
;  0 <= x < 1. From Abromowitz & Stegun. This is the integral of the function
;  f(x) = sqrt( 1 - (x*sin(t))^2 ) from t=0 to t=Pi/2. Fractional errors are
;  smaller than 10^(-8) when compared to Maple.
;
m = x^2
m1 = 1d - m
E = 1d + 0.44325141463d*m1     + 0.06260601220d*(m1^2) + $
         0.04757383546d*(m1^3) + 0.01736506451d*(m1^4) + $
       ( 0.24998368310d*m1     + 0.09200180037d*(m1^2) + $
         0.04069697526d*(m1^3) + 0.00526449639d*(m1^4) )*alog(1d/m1)

return, E
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lc_M0_S1half, beta
;
;  Calculate the m=0, s=1/2 laplace coefficient via an elliptic integral.
;  See equations (45-46), page 497, in Brouwer & Clemence (1961). Those
;  equations are only valid for beta<1, so the reciprocal relation
;  lc(beta) = alpha*lc(alpha) is used for alpha=1/beta when beta>1.
;
alpha = beta
j = where(beta gt 1, Nj)
if (Nj gt 0) then alpha[j] = 1/beta[j]
lc = (4/!dpi)*Elliptic_K(alpha)
if (Nj gt 0) then lc[j] = alpha[j]*lc[j]

return, lc
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lc_M1_S1half, beta
;
;  Calculate the m=1, s=1/2 laplace coefficient via elliptic integrals.
;  See equations (45-46), page 497, in Brouwer & Clemence (1961). Those
;  equations are only valid for beta<1, so the reciprocal relation
;  lc(beta) = alpha*lc(alpha) is used for alpha=1/beta when beta>1.
;
alpha = beta
j = where(beta gt 1, Nj)
if (Nj gt 0) then alpha[j] = 1/beta[j]
lc = (4/!dpi/alpha)*( Elliptic_K(alpha) - Elliptic_E(alpha) )
if (Nj gt 0) then lc[j] = alpha[j]*lc[j]

return, lc
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro lc_M0and1_s, s, beta, lc_M0_s, lc_M1_s
;
;  Calculate the s=1/2, 3/2, 5/2 (etc) laplace coefficients having m=0 (lc_M0_s)
;  and m=1 (lc_M1_s). This is done iteratively via equation (64) of Brouwer & 
;  Clemence (1961), page 501, with s replaced by s-1. Those equations are only
;  valid for beta<1, so the reciprocal relation lc(beta) = (alpha^(2s))*lc(alpha)
;  is used for alpha=1/beta when beta>1.


;Make argument < 1
alpha = beta
j = where(beta gt 1, Nj)
if (Nj gt 0) then alpha[j] = 1/beta[j]

;laplace coefficients with s=1/2 and m=0 and 1
lc_M0_Sonehalf = lc_M0_S1half(alpha)
lc_M1_Sonehalf = lc_M1_S1half(alpha)

;2*s = 1, 3, 5, etc
two_s = round(2*s)

;we are done if s=1/2
lc_M0_s = lc_M0_Sonehalf
lc_M1_s = lc_M1_Sonehalf

;Otherwise if s>1/2, calculate the laplace coefficients iteratively
if (two_s gt 1) then begin $

    sc = 1.5d & $
    two_sc = round(2*sc) & $
    lc_M0_Sminus1 = lc_M0_Sonehalf & $
    lc_M1_Sminus1 = lc_M1_Sonehalf & $
    while (two_sc le two_s) do begin $

        ss = two_sc/2d & $

        ;calculate laplace coefficient with m=0, s=ss+1
        m = 0 & $
        lc_M0_s = (ss - 1 - m)*(1 + alpha^2)*lc_M0_Sminus1 + $
                  2*(m + ss - 2)*alpha*lc_M1_Sminus1 & $
        factor = (ss - 1)*((1-alpha^2)^2) & $
        lc_M0_s = lc_M0_s/factor & $

        ;calculate laplace coefficient with m=1, s=ss+1
        m = 1 & $
        lc_M1_s = (ss - 1 - m)*(1 + alpha^2)*lc_M1_Sminus1 + $
                  2*(m + ss - 2)*alpha*lc_M0_Sminus1 & $
        factor = (ss - 1)*((1-alpha^2)^2) & $
        lc_M1_s = lc_M1_s/factor & $

        ;update lc_M0_Sminus1, lc_M1_Sminus1, and two_sc
        lc_M0_Sminus1 = lc_M0_s & $
        lc_M1_Sminus1 = lc_M1_s & $
        two_sc = two_sc + 2 & $

    endwhile & $
endif

;Use reciprocal relation if beta>1
if (Nj gt 0) then begin $
    lc_M0_s[j] = (alpha[j]^two_s)*lc_M0_s[j] & $
    lc_M1_s[j] = (alpha[j]^two_s)*lc_M1_s[j] & $
endif

return
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function lap_coeff, m, s, beta
;
;  Calculate the laplace coefficient iteratively via
;  equation (47) of Brouwer & Clemence (1961), page 497. Those
;  equations are only valid for beta<1, so the reciprocal relation
;  lc(beta) = (alpha^(2s))*lc(alpha) is used for alpha=1/beta when beta>1.
;  The laplace coefficient is the integral of
;  f(t)=(2/Pi)*cos(m*t)/(1+beta^2-2*beta*cos(t))^s integrated from t=0 to t=Pi,
;  where m=any integer, s=0.5, 1.5, 2.5, etc, and beta>0.0.

;Make argument < 1
alpha = beta
j = where(beta gt 1, Nj)
if (Nj gt 0) then alpha[j] = 1/beta[j]

;Calculate the m=0 and m=1 laplace coefficients
lc_M0and1_s, s, alpha, lc_M0, lc_M1

;Set laplace coefficient in case m=0 or 1
m_abs = abs(m)
if (m_abs eq 0) then lc_M = lc_M0
if (m_abs eq 1) then lc_M = lc_M1

;Use equation (47) to calculate the m>1 laplace coefficient iteratively. 
mc = 2
lc_Mminus2 = lc_M0
lc_Mminus1 = lc_M1
while (mc le m_abs) do begin $
    lc_M = (mc - 1)*(alpha + 1/alpha)*lc_Mminus1 - (mc + s - 2)*lc_Mminus2 & $
    lc_M = lc_M/(mc - s) & $
    lc_Mminus2 = lc_Mminus1 & $
    lc_Mminus1 = lc_M & $
    mc = mc + 1 & $
endwhile

;Use reciprocal relation if beta>1
if (Nj gt 0) then lc_M[j] = (alpha[j]^(2*s))*lc_M[j]

return, lc_M
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function deriv_lap_coeff, m, s, beta
;
;  Calculate the derivative of the laplace coefficient with respect to beta,
;  where m=integer, s=0.5, 1.5, 2.5, etc, and beta>0.0. This uses equation (67)
;  of Brouwer & Clemence (1961), page 502.

dlc = s*( lap_coeff(m-1, s+1, beta) - 2*beta*lap_coeff(m, s+1, beta) + $
          lap_coeff(m+1, s+1, beta) )

return, dlc
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

