-module(test).
-export [total/2, prefix_sum/3].
-import(wtree, [reduce/2, scan/5]).
-import(workers, [get/2, get/3, put/3]).

total(W, Key) ->
    N = 4,
    wtree:reduce(W,
        fun(PS) -> io:format("~p  ~n",[N]), lists:sum(workers:get(PS, Key)) end, %Leaf
        fun(L, R) ->  L + R end % combine
    ).

prefix_sum(W, Src, Dst) ->
    wtree:scan(W,
        fun(PS) -> lists:sum(workers:get(PS, Src)) end, %Leaf1 (going up)
        fun(PS, AccIn) ->
            workers:put(PS, Dst,
                element(1,
                    lists:mapfoldl(fun(E, Acc) -> {E+Acc, E+Acc} end,
                                    AccIn,
                                    workers:get(PS, Src)
                    )
                )
            ) end, % Leaf2 (going down)
        fun(L, R) -> L + R end, % combine
        0
    ).

%% Example:
% c(wtree).
% c(workers).
% c(test).
% W = wtree:create(4).
% workers:rlist(W, 15, 99, data).
% workers:retrieve(W, data).
% test:prefix_sum(W, data, data_acc).
% lists:append(workers:retrieve(W, data_acc)).
% lists:append(workers:retrieve(W, data)).


