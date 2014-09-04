srcs=$(wildcard *.cpp)
objs=$(srcs:%.cpp=%.o)
all: $(objs)
	g++ -o is_score.exe $(objs)
%.o : %.cpp
	echo $@
	g++ -o $@ -c $<
