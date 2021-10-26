
# https://www.rosettacode.org/wiki/Roman_numerals/Encode#Julia
function roman_encode(n::Integer)
    if n < 1 || n > 4999 throw(DomainError(n)) end
 
    DR = [["I", "X", "C", "M"] ["V", "L", "D", "MMM"]]
    rnum = ""
    for (omag, d) in enumerate(digits(n))
        if d == 0
            omr = ""
        elseif d <  4
            omr = DR[omag, 1] ^ d
        elseif d == 4
            omr = DR[omag, 1] * DR[omag, 2]
        elseif d == 5
            omr = DR[omag, 2]
        elseif d <  9
            omr = DR[omag, 2] * DR[omag, 1] ^ (d - 5)
        else
            omr = DR[omag, 1] * DR[omag + 1, 1]
        end
        rnum = omr * rnum
    end
    return rnum
end