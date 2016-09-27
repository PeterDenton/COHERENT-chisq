CXX=g++
CFlags=-c -Wall -O3 -pedantic -MMD -Werror
Sources=$(wildcard src/*.cpp)
IncludeDir=-I./include
AllObjects=$(addprefix obj/,$(notdir $(Sources:.cpp=.o)))
Executables=main
Objects=$(filter-out $(addprefix obj/,$(Executables:=.o)),$(AllObjects))

all: $(Sources) $(Executables)

$(Executables): $(AllObjects)
	@mkdir -p data obj
	$(CXX) $(CoverageFlags) $(Objects) $(addprefix obj/,$@.o) -o $@

obj/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFlags) $(CoverageFlags) $(IncludeDir) $< -o $@

-include $(AllObjects:.o=.d)

test: $(Executables)
	$(foreach exe,$(Executables),./$(exe);)

clean:
	rm -rf obj/*.o obj/*.d $(Executables)
