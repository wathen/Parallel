-module(path).

-export([path/0]).

% If you copy this to your own computer, set the directoy in the add_path
% to the location where you put the course library files.
path() -> code:add_path("/Users/michaelwathen/Dropbox/Parallel/Code/erl").
