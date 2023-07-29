#pragma once
int mapReduce(std::function<int(int, int)>  mapfunc, std::function<int(int, int)> redfunc, 
    const std::vector<int>& arg1,
    const std::vector<int>& arg2);

/*
template<class MAPFUNC, class REDFUNC>
auto mapReduce(MAPFUNC mapfunc, REDFUNC redfunc,
    
    */