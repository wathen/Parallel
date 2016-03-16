-module(hw2).

-export [pi/0, degree_to_radian/1, radian_to_degree/1, move_pos/2, move_par/4, q1/0].
-export [rle/1, rlep/2, longest_run/3, rleHasKey/2, longrun_combine/2, q3/0,  rlePos/2, rle_helper/2, listSplit/2, longrun_root/1].
-export [match_count/2, best_match/2, best_match/3, best_match_par/3, firstOccur/2, q4/0, best_match_combine/2, lrmatch/2, q4timing/0].
-export [q1_leaf1/2, combine/2].

% pi from Wolfram Alpha
pi() -> 3.1415926535897932384626433832795028841971693993751058.

degree_to_radian(Deg) -> Deg*pi()/180.0.

radian_to_degree(Rad) -> Rad*180/pi().


q1(W) ->
  InitialPos = {0.0,0.0,0},
  % Move = [{40,4}, {30,5}, {10,1}, {60,4}],
  Move = lists:map(fun(A) -> {10*A,A} end,[random:uniform(10) || _ <- lists:seq(1, 10)]),
  workers:update(W, data, misc:cut(Move, 2)),
  hw2:move_par(W, InitialPos, data, data_acc),
  io:format('~p~n',[lists:append(workers:retrieve(W, data_acc))]),
  B = move_pos(InitialPos, Move),
  io:format('~p~n',[B]).

  % L = lists:seq(1,10).


q1() ->
  W = wtree:create(2),
  q1(W).

q1_leaf1(ProcState, MoveKey) ->
    Src = workers:get(ProcState, MoveKey),
    move_pos({0.0,0.0,0.0},Src).

q1_leaf2(InitPos, Moves) ->
  element(1,
      lists:mapfoldl(
        fun(Move,Acc) -> {move_pos(Acc,Move), move_pos(Acc,Move)} end,
        InitPos,
        Moves
      )).
q1_combine(Left, Right) -> combine(Left, Right).

% move_par: parallel version of move_pos.
%   W is a worker-tree, and MoveKey is a key for ProcState.
%   The list of moves should be distributed across the workers of W and
%   associated with MoveKey.  The traveler starts at the position
%   indicated by InitPos.  The position after each move is stored
%   in the list associated with PosKey.  move_par returns the final
%   position of the traveller.
% move_par(W, InitPos, MoveKey, PosKey) -> % stub
%   use([W, InitPos, MoveKey, PosKey]),
%   X_final = 0.0,   % x-coordinate at final position
%   Y_final = 0.0,   % y-coordinate at final position
%   Dir_final = 0.0, % direction of travel at final direction
%                    %   (degrees counter-clockwise from the x-axis)
%   {X_final, Y_final, Dir_final}.move_par(W, InitPos, MoveKey, PosKey) -> % stub
move_par(W, InitPos, MoveKey, PosKey) ->
  Leaf1 = fun(ProcState) ->
    q1_leaf1(ProcState, MoveKey)
  end,
  Leaf2 = fun(PS, AccIn) ->
    workers:put(PS, PosKey, q1_leaf2(AccIn, workers:get(PS, MoveKey)))
  end,

  Combine = fun(Left, Right) ->
    q1_combine(Left, Right)
  end,
  wtree:scan(W, Leaf1, Leaf2, Combine, InitPos).

% move pos combine function
combine({X1,Y1,Dir1}, {X2,Y2,Dir2}) ->
  NewDir = Dir1+Dir2,
  Dist =  math:sqrt(math:pow(X2,2) + math:pow(Y2,2)),
  NewDirRad = math:atan2(Y2,X2)+degree_to_radian(Dir1),
  { X1 + Dist*math:cos(NewDirRad),
    Y1 + Dist*math:sin(NewDirRad),
    NewDir
  }.

% move_pos(InitPos, Move):
%  InitPos = {X, Y, Dir} where X and Y are the coordinates of the initial
%                          position.  Dir is the direction that the
%                          traveler is facing (in degrees counter-clockwise
%                          from the x-axis).
%  Move = {Dist, Dir}:  The traveler turns Dir degrees counter-clockwise,
%                          and then move Dist units "forward".
% We return the final position.
% If Move is a list, then we perform the given moves, in head-to-tail order.
move_pos({X, Y, Dir}, {Turn, Dist}) ->
  NewDir = Dir+Turn,
  NewDirRad = degree_to_radian(NewDir),
  { X + Dist*math:cos(NewDirRad),
    Y + Dist*math:sin(NewDirRad),
    NewDir
  };
move_pos(Pos, []) -> Pos;
move_pos(Pos, [MoveH | MoveTail]) ->
  move_pos(move_pos(Pos, MoveH), MoveTail).

% rle: run-length encoding
%   Convert a list of values in to a list of {Value, Count} pairs.
%     such that Count consecutive occurrences of Value are replaced by
%     the tuple {Value, Count}.  For example,
%     rle([1,2,2,3,3,3,0,5,5]) -> [{1,1},{2,2},{3,3},{0,1},{5,2}]
rle(L) when is_list(L) -> % stub
  rle_helper(L, []).

% Sequential longest run
seq_rle(L, Key) ->
  R = rlePos(rle_helper(L, []),0),
  K = lists:keysort(1,rleHasKey(R, Key)),
  if
    K =:= [] -> [];
    true -> lists:last(K)
  end.

% longest run - root function
rlep(L, Key) when is_list(L) -> % stub
  K = rlePos(rle_helper(L, []),0),
  listSplit(K, Key).

% longest run - root function helpers
rle_helper([ ] , Acc) -> lists:reverse(Acc);
rle_helper([H|T], []) -> rle_helper(T, [{H, 1}]);
rle_helper([H|T], [{Char,Count}|AT]) ->
    if H =:= Char -> rle_helper(T, [{Char, Count+1} | AT]);
        true -> rle_helper(T, [{H, 1} | [{Char, Count} | AT]])
    end.

rlePos([ ],_)->[];
rlePos([{Char,Count}| AT],0) -> [ {Char,Count,1} | rlePos(AT,1+Count)];
rlePos([{Char,Count}| AT],N) -> [ {Char,Count,N} | rlePos(AT,Count+N)].

rleHasKey(Tuple, Key) when is_tuple(Tuple) -> rleHasKey([Tuple], Key);
rleHasKey(List,Key) when is_list(List) ->
  [{Count, N} || {Char, Count, N} <- List, Char =:= Key].

listSplit([], _) -> [ 0, [], [], []];
listSplit([A], Key) ->
  [1,  rleHasKey(A, Key),  [], []];
listSplit([A,B],Key) ->
  H1 = rleHasKey(A, Key),
  T1 = rleHasKey(B, Key),
  [2,  H1,  T1, []];
listSplit(L=[H|M0],Key)  ->
  {M, T} = lists:split(length(M0) - 1, M0),
  H1 = rleHasKey(H, Key),
  T1 = rleHasKey(T, Key),
  M1 = rleHasKey(M, Key),
  [length(L)+2, H1, T1, M1].

% longest run - combine function
longrun_combine([N1, Ll, Lr, Lm], [N2, Rl, Rr, Rm]) ->
  lists:map(fun({C,A,B}) -> {C,A,B+N2} end, Rm),
  Match = longrun_match(Lr, Rl),
  M = lists:append(lists:append(Lm,Match),Rm),
  [N1+N2, Ll, Rr, M].

longrun_match([],[{A2,B2}]) -> [{A2,B2}];
longrun_match([{A1,B1}],[]) -> [{A1,B1}];
longrun_match([],[]) -> [];
longrun_match([{A1,B1}],[{A2,_}]) -> [{A1+A2,B1}].


% longest run - root function
longrun_root([_,F,L,M]) ->
  FinalX = lists:append(lists:append(F,M),L),
  SortX = lists:keysort(1,FinalX),
  if SortX =:= [] -> [];
     true -> lists:last(SortX)
  end.

% return the a description of the longest run of V's in the distributed
% list associated with Key in worker-pool W.
longest_run(W, Key, V) -> % stub
  wtree:reduce(W,
    fun(ProcState) -> rlep(workers:get(ProcState, Key), V) end,  % Leaf
    fun(Left, Right) -> longrun_combine(Left, Right) end, % Combine
    fun(X) -> longrun_root(X) end % Root
   ).

q3() ->
  W = wtree:create(4),
  Key = 1,
  MyList = lists:append([1,1,1,1,1,1], lists:seq(0,1000000)),
  workers:update(W, data, misc:cut(MyList, 4)),
  workers:retrieve(W, data),

  ParTime = time_it:t(fun() -> longest_run(W, data, Key) end,20),
  Time = time_it:t(fun() -> seq_rle(MyList, Key) end,20),
  io:format("~p  ~n",[Time]),
  io:format("~p  ",[ParTime]).


timing(P, [H]) ->
  W = wtree:create(P),
  L2 = lists:seq(0,H),
  workers:update(W, data, misc:cut(L2, P)),
  L1 = lists:seq(400,600),
  ParTime = time_it:t(fun() -> best_match_par(W, L1, data) end,20),
  Time = time_it:t(fun() -> best_match(L1,L2) end,20),
  STime = element(2,hd(Time)),
  PTime = element(2,hd(ParTime)),
  io:format(" ~6b  ~6b  ~6b  ~10.4e  ~10.4e   ~10.4e        ~10.4e~n", [P, length(L1), H, STime, PTime, STime/(length(L1)*H), PTime*(P/(length(L1)*H))]);
timing(P, [H|T]) ->
  W = wtree:create(P),
  L2 = lists:seq(0,H),
  workers:update(W, data, misc:cut(L2, P)),
  L1 = lists:seq(400,600),
  ParTime = time_it:t(fun() -> best_match_par(W, L1, data) end,20),
  Time = time_it:t(fun() -> best_match(L1,L2) end,20),
  STime = element(2,hd(Time)),
  PTime = element(2,hd(ParTime)),
  io:format(" ~6b  ~6b  ~6b  ~10.4e  ~10.4e   ~10.4e        ~10.4e~n", [P, length(L1), H, STime, PTime, STime/(length(L1)*H), PTime*(P/(length(L1)*H))]),
  % io:format("Sequential ~p  ~n",[Time]),
  % io:format("Parallel   ~p  ~n",[ParTime]),
  timing(P,T).


q4timing() ->
  io:format("      P     N_1    N_2       T_s         T_p      T_s/(N_1*N_2)     T_p*P/(N_1*N_2)~n"),
  timing(4,[15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000, 100000]).




q4() ->
W = wtree:create(4),
MyList = lists:seq(0,10000000),
workers:update(W, data, misc:cut(MyList, 4)),
List1 = lists:seq(3,6),

best_match_par(W, List1, data).

% match_count:
%   We return the number of values for I,
%   with 1 <= I <= min(length(L1), length(L2)), such that
%   lists:nth(I, L1) =:= lists:nth(I, L2).
%   A brute force way to do this would be
%     length([ok || I <- lists:seq(1, min(length(L1), length(L2))),
%                        lists:nth(I, L1) =:= lists:nth(I, L2)])
%   but that would take quadratic time in the length of the lists :(
%   This implementation is linear time.
match_count(L1, L2) when is_list(L1), is_list(L2) -> match_count(L1, L2, 0).
match_count([], _, C) -> C;
match_count(_, [], C) -> C;
match_count([H | T1], [H | T2], C) -> match_count(T1, T2, C+1);
match_count([_ | T1], [_ | T2], C) -> match_count(T1, T2, C).

firstOccur([{H1,T1}],[_]) -> {H1,T1};
firstOccur([{H1,T1}],[]) -> {H1,T1};
firstOccur(L,[]) ->
        {H,T} = lists:split(length(L) - 1, L),
        firstOccur(H,T);
firstOccur(H,T) ->
        % io:format('H~p~n',[H]),
        {H1,T1} = lists:split(length(H) - 1, H),
        [{M1, A1}] = T,
        [{M2, _}] = T1,
        if
           M2 < M1 -> {M1, A1};
           M2 =:= M1, M1 > 0 -> firstOccur(H1,T1);
           true -> {M1,A1}
        end.

tupleAdd([{A,B}],N) -> [{A,B+N}];
tupleAdd([{A,B}|T], N) -> [{A,B+N}| tupleAdd(T,N)].

best_match_leaf(L1,L2) ->
  N = length(L1)-1,
  X = best_match(L1, L2, -N),
  LengthL1 = length(L1),
  LengthL2 = length(L2),
  if
    LengthL1 >= LengthL2 -> [{LengthL1,LengthL2}, X];
    LengthL1 < LengthL2 ->
          {H, M1} = lists:split(LengthL1-1, X),
          {M, T} = lists:split(length(M1)-LengthL1+1, M1),
          SortX = lists:keysort(1,M),
          Mfinal = firstOccur(SortX,[]),
          [{LengthL1,LengthL2}, H, T, Mfinal];
    true -> []
  end.

best_match_combine([{N,Ln}, Ll, Lr, {Lm1,Lm2}], [{N,Rn}, Rl, Rr, {Rm1,Rm2}]) ->
  MatchLR = best_match_LR_combine(Ln,Lr,Rl),
  M = lists:append(lists:append([{Lm1,Lm2}],MatchLR),[{Rm1,Rm2+Ln}]),
  Mmax = firstOccur(lists:keysort(1,M),[]),
  R = tupleAdd(Rr,Ln),
  [{N,Rn+Ln}, Ll, R, Mmax];
best_match_combine([{N,Ln}, L], [{N,Rn}, R]) ->
  MatchLR = best_match_LR_combine(Ln,L,R),
  if
    N < Ln+Rn ->
          {H, M1} = lists:split(N-1, MatchLR),
          {M, T} = lists:split(length(M1)-N+1, M1),
          SortX = lists:keysort(1,M),
          Mfinal = firstOccur(SortX,[]),
          [{N,Ln+Rn}, H, T, Mfinal];
    true -> [{N,Ln+Rn}, MatchLR]
  end.


best_match_LR_combine(N, L, R) -> lrmatch(L,tupleAdd(R,N)).

lrmatch(F,[]) -> F;
lrmatch([],F) -> F;
lrmatch([{M1,A1}|T1],[{M2,A1}|T2]) -> [{M1+M2,A1} | lrmatch(T1,T2)];
lrmatch([{M1,A1}|T1],[{M2,A2}|T2]) -> [{M1,A1} | lrmatch(T1,[{M2,A2}|T2])].

best_match_root([{_,_}, L, R, M]) ->
  Final = lists:append(lists:append(L,[M]),R),
  firstOccur(lists:keysort(1,Final),[]).


% best_match(L1, L2) -> {MatchCount, Alignment}
%   Find the smallest value of Alignment, with
%   -length(L1) =< Alignment =< length(L2) that maximizes
%     MatchCount = if
%       Alignment <  0 -> match_count(lists:nthtail(-Alignment, L1), L2);
%       Alignment >= 0 -> match_count(L1, lists:nthtail(L2))
%     end
%   Examples:
%     best_match([1,2,3],[1,0,2,3,2,1,2,0]) -> {2,1}
%     best_match("banana", "ananas is the genus of pineapple") -> {5, -1}
%     best_match("banana", "bandanna") -> {3, 0}
best_match(L1, L2) when is_list(L1), is_list(L2) ->
  N = length(L1)-1,
  SortX = lists:keysort(1,best_match(L1, L2, -N)),
  firstOccur(SortX,[]). % stub

best_match(_,[], _) -> [];
best_match(L1,L2, Alignment) ->
    if
      Alignment < 0 ->
                MatchCount = match_count(lists:nthtail(-Alignment, L1), L2),
                [{MatchCount, Alignment} | best_match(L1,L2, Alignment+1)];
      Alignment >= 0 ->
                MatchCount =  match_count(L1, L2),
                [{MatchCount, Alignment} | best_match(L1,tl(L2), Alignment+1)];
      true -> [best_match(L1,L2, Alignment+1)]
    end.


% best_match_par(W, Key1, Key2) -> {MatchCount, Alignment}
%   The parallel version of best_match.
%   best_match_par(W, Key1, Key2) should return the same value as
%   best_match(workers:retrieve(W, Key1), workers:retrieve(W, Key2))
%   but best_match should do it's work in parallel.
best_match_par(W, L1, L2) -> % stub
  wtree:reduce(W,
    fun(ProcState) -> best_match_leaf(L1, workers:get(ProcState, L2)) end,  % Leaf
    fun(Left, Right) -> best_match_combine(Left, Right) end, % Combine
    fun(X) -> best_match_root(X) end % Root
   ).


