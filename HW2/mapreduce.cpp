#include <iostream>
#include <tbb/tbb.h>

using namespace std;


int mapReduce(std::function<int(int, int)>  mapfunc, std::function<int(int, int)> redfunc, 
    const std::vector<int>& arg1,
    const std::vector<int>& arg2){    // Create an instance of tbb::combinable to hold the tmp results

    tbb::combinable<int> tmp_results;
    int dot = 0;

    tbb::parallel_for(tbb::blocked_range<size_t>(0, arg1.size()),
                        [&](const tbb::blocked_range<size_t>& r) {
                            for(size_t i = r.begin(); i!=r.end(); i++) {
                                dot += mapfunc(arg1[i], arg2[i]);
                            }
                            tmp_results.local() = dot;
                        }
    );

    return tmp_results.combine(redfunc);
}


int main (){
    std::vector<int> arg1 = {2, 3, 6, 1, 2, 0, 3, 2};
    std::vector<int> arg2 = {0, 4, 5, 1, 0, 2, 3, 6};

    // mapfunc
    std::function<int(int, int)> mapfunc = [](int x_i, int y_i) {
            return x_i * y_i;
    };

    // redfunc
    std::function<int(int, int)> redfunc = [](int x, int y) {
        return x + y;
    };

    int dot_product_result = mapReduce(mapfunc, redfunc, arg1, arg2);

    std::cout << "The dot product of the two vectors is - " << dot_product_result << std::endl;

    return 0;
}