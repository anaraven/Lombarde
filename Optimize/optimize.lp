#begin_lua
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

#const vshapesFile="lombarde.out".

five(X,Y,Z,W,A) :- (X,Y,Z,W,A) := @readFive(vshapesFile).

% Input: optimExplanation(I,A,B) when explanation I is in Optim(A,B)
optimExplanation(I,A,B) :- five("vshape",I,A,B,_).

% Input: arcInExplanation(I,X,Y,W) when explanation I has an arc (X,Y) with cost W 
arcInExplanation(I,X,Y,W) :- five("arcInVshape",I,X,Y,W).

% coexp(X,Y) means that {X,Y} are co-expressed.
% It is implicit in the definition of optimExplanation(_,X,Y).
coexp(X,Y) :- optimExplanation(_,X,Y).

% only one explanation is chosen for each {t,t'} in O 
1{ chosen(I) : optimExplanation(I,A,B) }1 :- coexp(A,B).

% consider the arcs that are part of the chosen explanation
chosenArc(X,Y,W) :- arcInExplanation(I,X,Y,W), chosen(I).

% minimize the global cost
#minimize [chosenArc(_,_,W) = W].
#hide.
#show chosenArc/3.

