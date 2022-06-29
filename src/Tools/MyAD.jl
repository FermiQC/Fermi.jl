import Base: +, -, *, ^, exp, sin

struct Dual <: Number
    x
    dx
end

+(a::Real, b::Dual) = Dual(a + b.x, b.dx)

+(a::Dual, b::Dual) = Dual(a.x + b.x, a.dx + b.dx)

*(a::Real, b::Dual) = Dual(a*b.x, a*b.dx)

*(a::Dual, b::Dual) = Dual(a.x * b.x, a.dx * b.x + a.x * b.dx)

exp(a::Dual) = Dual(exp(a.x), exp(a.x) * a.dx)

sin(a::Dual) = Dual(sin(a.x), cos(a.x) * a.dx)



#import Base: fill!, real, conj, float, abs, isless, exp
#
#
#function Dual(x::Real)
#    Dual(x,1.0)
#end
#
#function float(x::Dual)
#    return x
#end
#
#function conj(x::Dual)
#    return x
#end
#
#function real(x::Dual)
#    return x
#end
#
#function abs(a::Dual)
#    return Dual(abs(a.x), a.dx)
#end
#
#function isless(a::Dual, b::Real)
#    return isless(a.x, b)
#end
#
#function isless(b::Real, a::Dual)
#    return isless(a.x, b)
#end
#
#function isless(a::Dual, b::Dual)
#    return isless(a.x, b.x)
#end
#
#function exp(a::Dual)
#    return Dual(exp(a.x), exp(a.dx)*a.dx)
#end
#
#function *(a::Real, b::Dual)
#    return Dual(a*b.x, a*b.dx)
#end
#
#function *(a::Dual, b::Real)
#    return Dual(b*a.x, b*a.dx)
#end
#
#function *(a::Dual, b::Dual)
#    return Dual(a.x * b.x, a.dx * b.x + a.x * b.dx)
#end
#
#function +(a::Real, b::Dual)
#    return Dual(a + b.x, b.dx)
#end
#
#function +(b::Dual, a::Real)
#    return Dual(a + b.x, b.dx)
#end
#
#function +(a::Dual, b::Dual)
#    return Dual(a.x + b.dx, a.dx + b.dx)
#end
#
#function -(a::Real, b::Dual)
#    return Dual(a - b.x, b.dx)
#end
#
#function -(b::Dual, a::Real)
#    return Dual(a.x - b.x, b.dx)
#end
#
#function -(a::Dual, b::Dual)
#    return Dual(a.x - b.dx, a.dx - b.dx)
#end
#
#function fill!(dest::Matrix{Dual}, x::Float64)
#    dest .= Dual(x,1.0)
#end