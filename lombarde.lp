#begin_lua
function readCoexp(f)
  local file = io.open(f)
  local res = {}
  local pat = "(%S+)%s+(%S+)"
  while true do
    local x, y
    local row = file:read()
    if row == nil then break end
    for x, y in string.gfind(row, pat) do
      res[#res + 1] = Val.new(Val.FUNC, {x,y})
    end
  end
  file:close()
  return res
end

function readFive(f)
  local file = io.open(f)
  local res = {}
  local pat = "(%S+)%s+(%S+)%s+(%S+)%s+(%S+)%s+(%S+)"
  while true do
    local x, y, z, w, a
    local row = file:read()
    if row == nil then break end
    for x, y, z, w, a in string.gfind(row, pat) do
      res[#res + 1] = Val.new(Val.FUNC, {x,y,z,w, Val.new(Val.NUM,a)})
    end
  end
  file:close()
  return res
end
#end_lua.

#const vshapesFile="Out1/RDB8.1/mi-pearson-c3-133/graph-operons-3.out".

five(X,Y,Z,W,A) :- (X,Y,Z,W,A) := @readFive(vshapesFile).
vShape(I,X,Y) :- five("vshape",I,X,Y,_).
arcInVshape(I,X,Y,W) :- five("arcInVshape",I,X,Y,W).
coexp(X,Y) :- vShape(_,X,Y).

1{ chosen(I) :  vShape(I,A,B) }1 :- coexp(A,B).

used_arc(X,Y,W) :- arcInVshape(I,X,Y,W), chosen(I).

#minimize [used_arc(_,_,W) = W].
#hide.
#show used_arc/3.
