#begin_lua
function readThree(f)
  local file = io.open(f)
  local res = {}
  local pat = "(%S+)%s+(%S+)%s+(%S+)"
  while true do
    local x, y
    local row = file:read()
    if row == nil then break end
    for x, y, a in string.gfind(row, pat) do
      res[#res + 1] = Val.new(Val.FUNC, {x, y, Val.new(Val.NUM, a)})
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

#const firstFile="Out1/RDB8.1/M3D/pearson/mrnet/100K/Prodoric/MEME/gl-7.ncol".
#const secondFile="Out1/RDB8.1/M3D/pearson/mrnet/100K/Prodoric/MEME-GS-0/gl-7.ncol".
#const thirdFile="Input/RDB8.1/gold-std.ncol".

arc(net1,X,Y,W) :- (X,Y,W) := @readThree(firstFile).
arc(net2,X,Y,W) :- (X,Y,W) := @readThree(secondFile).
arc(net3,X,Y,W) :- (X,Y,W) := @readThree(thirdFile).

arc(N,X,Y) :- arc(N,X,Y,W).
intersect(A,B) :- arc(net1,A,B), arc(net2,A,B), arc(net3,A,B).
only_in(net1,A,B) :- arc(net1,A,B), not arc(net2,A,B), not arc(net3,A,B).
only_in(net2,A,B) :- arc(net2,A,B), not arc(net1,A,B), not arc(net3,A,B).
only_in(net3,A,B) :- arc(net3,A,B), not arc(net1,A,B), not arc(net2,A,B).
both(net1,net2,A,B) :- arc(net1,A,B), arc(net2,A,B), not arc(net3,A,B).
both(net1,net3,A,B) :- arc(net1,A,B), arc(net3,A,B), not arc(net2,A,B).
both(net2,net3,A,B) :- arc(net2,A,B), arc(net3,A,B), not arc(net1,A,B).

count(A,B,C,D,E,F,G) :- A=[only_in(net1,_,_)], B=[only_in(net2,_,_)], C=[only_in(net3,_,_)], 
	D=[both(net1,net2,_,_)], E=[both(net1,net3,_,_)], F=[both(net2,net3,_,_)], 
	G=[intersect(_,_)].

r(0,0,0,0).
r(0,0,1,A) :- A=[only_in(net1,_,_)].
r(0,1,0,B) :- B=[only_in(net2,_,_)].
r(1,0,0,C) :- C=[only_in(net3,_,_)].
r(0,1,1,D) :- D=[both(net1,net2,_,_)].
r(1,0,1,E) :- E=[both(net1,net3,_,_)].
r(1,1,0,F) :- F=[both(net2,net3,_,_)].
r(1,1,1,G) :- G=[intersect(_,_)].
#hide.
% #show intersect/2.
% #show only_in/3.
% #show both/4.
#show count/7.
#show r/4.
