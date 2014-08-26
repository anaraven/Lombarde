#begin_lua
function readTriplet(f)
  local file = io.open(f)
  local res = {}
  local pat = "(%S+)%s+(%S+)%s+(%S+)"
  while true do
    local x, y, z
    local row = file:read()
    if row == nil then break end
    for x, y, z in string.gfind(row, pat) do
      res[#res + 1] = Val.new(Val.FUNC, {x,y, Val.new(Val.NUM,z)})
    end
  end
  file:close()
  return res
end
#end_lua.

#const coexpFile="Input/RDB8.1/M3D/pearson/mrnet/small.ncol".
#const arcsFile="Input/RDB8.1/Prodoric/MEME/small.ncol".

coexp(X,Y) :- (X,Y,W) := @readTriplet(coexpFile).
arc(X,Y,10**W) :- (X,Y,W) := @readTriplet(arcsFile).

{ used_arc(X,Y,W) } :- arc(X,Y,W).
precedes(X,Y) :- used_arc(X,Y,_).
precedes(X,Y) :- precedes(X,Z), used_arc(Z,Y,_).

v_connected(A,B) :- precedes(M,A), precedes(M,B), coexp(A,B).

:- coexp(A,B), not v_connected(A,B).

#minimize [used_arc(X,Y,W)=W].

#hide.
#show used_arc/3.
