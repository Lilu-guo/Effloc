source=glocate27.cpp divsufsort.o sssort.o trsort.o loadkit.o savekit.o InArray.cpp #加速，减少new，重复使用
# source=lookBit.cpp                               #查看索引中的位
# source=rpRead2Pattern.cpp                        #将rpRead转成Original_s的patterns格式

all: clean effloc
effloc:$(source)
# samStatics:$(source)
	# g++ -pg -std=c++11 -o guo $(source) #gprofile
	# g++ -g -std=c++11 -o guo $(source)  #gdb
	# g++ -std=c++11 -o guo $(source)
	# g++ -DMYDEBUG -g -std=c++11 -o $@ $(source) -lpthread                                        #gdb
	# g++ -DMYDEBUG -g -std=c++11 -L ./divsufsort/ -o $@ $(source) -ldivsufsort -lpthread          #gdb
	g++ -DMYDEBUG -O3 -L ./divsufsort/ -o $@ $(source) -ldivsufsort -mpopcnt -Winline     #加-O3优化选项
clean:
	rm -rf effloc