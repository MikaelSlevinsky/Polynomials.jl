immutable RatPoly{T<:Number}
    p::Poly{T}
    q::Poly{T}
    var::Symbol
    function RatPoly(p::Poly{T},q::Poly{T})
        @assert_variable(p,q)
        d = gcd(p, q)
        pd = div(p, d)
        qd = div(q, d)
        n = findfirst(qd)
        qdn = qd[n]
        if abs(qdn) ≤ 2eps(T)
            new(pd, qd, p.var)
        else
            new(pd/qdn, qd/qdn, p.var)
        end
    end
end
RatPoly{T<:Number}(p::Poly{T}, q::Poly{T}) = RatPoly{T}(p,q)
RatPoly{T<:Integer}(p::Poly{T}, q::Poly{T}) = RatPoly{Rational{T}}(convert(Poly{Rational{T}},p),convert(Poly{Rational{T}},q))

eltype{T}(::RatPoly{T}) = T

RatPoly{T<:Number,S<:Number}(p::Poly{T}, q::Poly{S}) = RatPoly(promote(p,q)...)


convert{S,T<:Number}(::Type{RatPoly{S}},P::T) = convert(Poly{S},P)/convert(Poly{S},one(P))
convert{S,T}(::Type{RatPoly{S}},P::Poly{T}) = convert(Poly{S},P)/convert(Poly{S},one(P))
convert{S<:Integer,T<:Integer}(::Type{RatPoly{S}},P::T) = convert(Poly{S},P)//convert(Poly{S},one(P))
convert{S<:Integer,T<:Integer}(::Type{RatPoly{S}},P::Poly{T}) = convert(Poly{S},P)//convert(Poly{S},one(P))

promote_rule{S,T<:Number}(::Type{RatPoly{S}},::Type{T}) = RatPoly{promote_type(S,T)}
promote_rule{S,T}(::Type{RatPoly{S}},::Type{Poly{T}}) = RatPoly{promote_type(S,T)}


ratpolyval(r::RatPoly,x) = polyval(r.p,x)./polyval(r.q,x)
ratpolyval{T<:Union{Rational,Integer},S<:Union{Rational,Integer}}(r::RatPoly{T},x::S) = polyval(r.p,x).//polyval(r.q,x)

if VERSION >= v"0.4"
    call(r::RatPoly, x) = ratpolyval(r, x)
end

copy(r::RatPoly) = RatPoly(copy(r.p), copy(r.q))

zero(r::RatPoly) = zero(r.p)//one(r.q)
zero{T}(::Type{RatPoly{T}}) = zero(Poly{T})//one(Poly{T})
one(r::RatPoly) = one(r.p)//one(r.q)
one{T}(::Type{RatPoly{T}}) = one(Poly{T})//one(Poly{T})


//(p::Poly, q::Poly) = RatPoly(p, q)
//(R1::RatPoly, R2::RatPoly) = R1/R2

+(R::RatPoly) = RatPoly(+R.p,R.q)
-(R::RatPoly) = RatPoly(-R.p,R.q)
+(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q+R2.p*R1.q,R1.q*R2.q)
-(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q-R2.p*R1.q,R1.q*R2.q)
*(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.p,R1.q*R2.q)
/(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q,R1.q*R2.p)

for op in (:+,:-,:*,:/,://)
    @eval begin
        $op(R::RatPoly, P::Union{Poly,Number}) = $op(promote(R,P)...)
        $op(P::Union{Poly,Number}, R::RatPoly) = $op(promote(P,R)...)
    end
end

# compute the (m,n)-Pade approximant to the Polynomial with coefficients c.

function pade{T}(c::Poly{T},m::Int,n::Int)
    @assert m+n < length(c)
    rold,rnew = Poly([zeros(T,m+n+1);one(T)],c.var),Poly(c.a[1:m+n+1],c.var)
    uold,vold = Poly([one(T)],c.var),Poly([zero(T)],c.var)
    unew,vnew = vold,uold
    for i=1:n
        temp0,temp1,temp2 = rnew,unew,vnew
        q,rnew = divrem(rold,rnew)
        unew,vnew = uold-q*unew,vold-q*vnew
        rold,uold,vold = temp0,temp1,temp2

    end
    if vnew[0] == 0
        d = gcd(rnew,vnew)
        rnew /= d
        vnew /= d
    end
    RatPoly(rnew/vnew[0],vnew/vnew[0])
end
