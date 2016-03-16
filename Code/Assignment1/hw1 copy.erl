-module(hw1).
-export([nthtail/2,prefix/2,search/2,
         subPrint/2,subtract/2]).


-spec nthtail(N, List) -> Tail when
      N :: non_neg_integer(),
      List :: [T],
      Tail :: [T],
      T :: term().

-spec prefix(List1, List2) -> boolean() when
      List1 :: [T],
      List2 :: [T],
      T :: term().

-spec search(List1, List2, N) -> List3 when
      N :: non_neg_integer(),
      List1 :: [T],
      List2 :: [T],
      List3 :: [T],
      T :: term().

-spec search(List1, List2) -> List3 when
      List1 :: [T],
      List2 :: [T],
      List3 :: [T],
      T :: term().

-spec subtract(List1, List2) -> List3 when
      List1 :: [T],
      List2 :: [T],
      List3 :: [T],
      T :: term().

% EXPORTED FUNCTION
% nthtail function: question 1a)
% inputs :
%         N -> integer > 0
%         T -> List
% outputs :
%         R -> List
%
% Example usage:
% hw1:nthtail(1, [1,2,3]).
% Output: [2,3]
nthtail(1, [_|T]) -> T;
nthtail(0, L) -> L;
nthtail(N, [_|T]) -> nthtail(N - 1, T).

% EXPORTED FUNCTION
% prefix function: question 1b)
% inputs :
%         Head -> List
%         Tail -> List
% outputs :
%         R -> boolean
%
% Example usage:
% hw1:prefix([2], [1,2,3]).
% Output: false
prefix([X|Head], [X|Tail]) -> prefix(Head, Tail);
prefix([], List) -> true;
prefix([_|_], List) -> false.

% EXPORTED FUNCTION
% search/2 function: question 1c)
% inputs :
%         L1 -> List
%         L2 -> List
% outputs :
%         R -> List
%
% Example usage:
% hw1:search([2], [1,2,3]).
% Output: [2]
search(L1, L2) -> search(L1, L2, 1).

% HELPER FUNCTION
% search/3 function: helper function for question 1c)
% inputs :
%         L1 -> List
%         L1 -> List
%         Pos -> interger
% outputs :
%         R -> List
search(L1, L2, Pos) ->
    case {prefix(L1,L2), L2} of
         {true, [_|T]} -> [Pos | search(L1, T, Pos+1)];
         {true, []} -> [Pos];
         {false, []} -> [];
         {false, [_|T]} -> search(L1, T, Pos+1)
     end.

% EXPORTED FUNCTION
% subtract/2 function: question 2c)
% inputs :
%         L1 -> List
%         L1 -> List
% outputs :
%         R -> List
%
% Example usage:
% hw1:subtract([1,2], [1,2,3]).
% Output: []
subtract([],List) -> [];
subtract(List,[]) -> List;
subtract(L1,L1) -> [];
subtract(L1,L2) when is_list(L1),is_list(L2) ->
    SortedL1 = lists:keysort(1,lists:zip(L1,lists:seq(1,length(L1)))),
    SortedL2 = lists:sort(L2),
    F = elementCheck(SortedL1,SortedL2),
    {Result,_} = lists:unzip(lists:keysort(2,F)),
    Result.


% HELPER FUNCTION
% elementCheck/2 function: question 2c)
% inputs :
%         L1 -> List
%         L1 -> List
% outputs :
%         R -> List
elementCheck([],List) -> [];
elementCheck(List,[]) -> List;
elementCheck(H, [H2 | T2]) ->
    [{P1,P2}|T1] = H,

    if P1 > H2 -> [{P1,P2} | elementCheck(H, [T2])];
       P1 < H2 -> [{P1,P2} | elementCheck(T1, [H2|T2])];
       P1=:=H2 -> elementCheck(T1, T2)
    end.


% HELPER FUNCTION
% libSubtraction/1 function: question 2a)
% inputs :
%         N -> interger
% outputs :
%         R -> tuple
libSubtraction(N)->
    L1 = lists:seq(1,N),
    L2 = lists:reverse(L1),
    timer:tc(lists,subtract,[L1,L2]).


% HELPER FUNCTION
% customSubtraction/1 function: used in printing for 2d)
% inputs :
%         N -> interger
% outputs :
%         R -> tuple
customSubtraction(N)->
    L1 = lists:seq(1,N),
    L2 = lists:reverse(L1),
    timer:tc(hw1,subtract,[L1,L2]).

% EXPORTED FUNCTION
% subPrint/2 function: printing function 2d)
% inputs :
%         N -> interger
%         N -> list
%
% example usage:
% D = [500,100].
% subPrint(length(D),D).
% Output:
% 500  8  0.0032  58  0.2965217736169845  1.16
% 100  46  0.0046  18  0.039086503371292665  0.18
subPrint(0,[ ])->[];
subPrint(N,[H|T]) ->
    F1 = libSubtraction(H),
    [W1,[]] = tuple_to_list(F1),
    io:format("~p  ",[H]),
    io:format("~p  ",[W1]),
    io:format("~p  ",[W1/(H*H)]),
    F2 = customSubtraction(H),
    [W2,[]] = tuple_to_list(F2),
    io:format("~p  ",[W2]),
    io:format("~p  ",[W2/(H*math:log(H))]),
    io:format("~p  ~n",[W2/H]),
    subPrint(N-1,T).
