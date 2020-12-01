CFLAGS    :=`root-config --cflags --libs`

DllSuf    := so
SrcSuf    := c
ObjSuf    := o

INCLUDES  := -I./

CXXFLAGS  += $(INCLUDES) -std=c++11 -fPIC -O3 -Wall -Wpedantic

CLASSIFICATION_TOOLS_LIB := libclassification_tools.$(DllSuf)
SRCS = $(wildcard *.$(SrcSuf))
OBJS = $(patsubst %.$(SrcSuf), %.$(ObjSuf), $(SRCS))

.$(SrcSuf).$(ObjSuf):
	@$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -c $< $(SYSLIB) $(CFLAGS)"

all: $(CLASSIFICATION_TOOLS_LIB)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

$(CLASSIFICATION_TOOLS_LIB): $(OBJS)
	@$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS)
	@echo "$(CXX) $(CXXFLAGS) -shared -o $@ $^ $(CFLAGS)"

install:

.PHONY: distclean
distclean:
	@rm -f $(CLASSIFICATION_TOOLS_LIB) $(OBJS)

.PHONY: clean
clean:
	@$(RM) -f $(CLASSIFICATION_TOOLS_LIB) $(OBJS)

.PHONY: lint
lint:
	$(LINT) $(INC_SRCH_PATH) $(SRCS)
