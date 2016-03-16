-module(test).          % module attribute
-export([b/1]).   % module attribute

b(N)->
    L1 = lists:seq(1,N),
    L2 = lists:reverse(L1),
    timer:tc(lists,subtract,[L1,L2]).