immutable RatPoly{T<:Number}
    p::Poly{T}
    q::Poly{T}
    var::Symbol
    function RatPoly(p::Poly{T},q::Poly{T})
        @assert_variable(p,q)
        d = gcd(p, q)
        pd = div(p, d)
        qd = div(q, d)
        qd0 = qd[0]
        if abs(qd0) ≤ eps(T)
            new(pd, qd, p.var)
        else
            new(pd/qd0, qd/qd0, p.var)
        end
    end
end
RatPoly{T<:Number}(p::Poly{T}, q::Poly{T}) = RatPoly{T}(p,q)

eltype{T}(::RatPoly{T}) = T

RatPoly{T<:Number,S<:Number}(p::Poly{T}, q::Poly{S}) = RatPoly(promote(p,q)...)


Base.convert{S,T<:Number}(::Type{RatPoly{S}},P::T) = convert(Poly{S},P)/convert(Poly{S},one(P))
Base.convert{S,T}(::Type{RatPoly{S}},P::Poly{T}) = convert(Poly{S},P)/convert(Poly{S},one(P))

Base.promote_rule{S,T<:Number}(::Type{RatPoly{S}},::Type{T}) = RatPoly{promote_type(S,T)}
Base.promote_rule{S,T}(::Type{RatPoly{S}},::Type{Poly{T}}) = RatPoly{promote_type(S,T)}


rateval(r::RatPoly,x) = polyval(r.p,x)./polyval(r.q,x)

Base.call(r::RatPoly,x) = rateval(r,x)

copy(r::RatPoly) = RatPoly(copy(r.p), copy(r.q))

zero(r::RatPoly) = zero(r.p)/one(r.q)
zero{T}(::Type{RatPoly{T}}) = zero(Poly{T})/one(Poly{T})
one(r::RatPoly) = one(r.p)/one(r.q)
one{T}(::Type{RatPoly{T}}) = one(Poly{T})/one(Poly{T})


/(p::Poly, q::Poly) = RatPoly(p, q)

+(R::RatPoly) = RatPoly(+R.p,R.q)
-(R::RatPoly) = RatPoly(-R.p,R.q)
+(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q+R2.p*R1.q,R1.q*R2.q)
-(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q-R2.p*R1.q,R1.q*R2.q)
*(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.p,R1.q*R2.q)
/(R1::RatPoly,R2::RatPoly) = RatPoly(R1.p*R2.q,R1.q*R2.p)

+(R::RatPoly,P::Union{Poly,Number}) = +(promote(R,P)...)
+(P::Union{Poly,Number},R::RatPoly) = +(promote(P,R)...)
-(R::RatPoly,P::Union{Poly,Number}) = -(promote(R,P)...)
-(P::Union{Poly,Number},R::RatPoly) = -(promote(P,R)...)

*(R::RatPoly,P::Union{Poly,Number}) = *(promote(R,P)...)
*(P::Union{Poly,Number},R::RatPoly) = *(promote(P,R)...)
/(R::RatPoly,P::Union{Poly,Number}) = /(promote(R,P)...)
/(P::Union{Poly,Number},R::RatPoly) = /(promote(P,R)...)



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
