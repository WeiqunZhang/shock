
DEBUG = FALSE
PRECISION = DOUBLE
FSANITIZE = FALSE

CXX ?= g++
DEPFLAGS = -MMD -MP
WARNFLAGS = -Wall -Wextra -Wpedantic

ifeq ($(DEBUG),TRUE)
  BUILD_TYPE = .debug
  CXXFLAGS ?= -g -O0
  CPPFLAGS += -DDEBUG
else
  BUILD_TYPE = .opt
  CXXFLAGS ?= -g1 -O3
endif

CXXFLAGS += $(WARNFLAGS)

# Detect compiler type and set flags accordingly
ifneq (,$(shell $(CXX) --version 2>/dev/null | grep -i clang))
  CXX_TYPE = .clang
  # CXXFLAGS += -stdlib=libc++    # use libc++ if libstdc++ is not found
else ifneq (,$(shell $(CXX) --version 2>/dev/null | grep -i 'Free Software Foundation'))
  CXX_TYPE = .gcc
else
  CXX_TYPE =
endif

ifeq ($(PRECISION),FLOAT)
  PREC_TYPE = .float
  CPPFLAGS += -DUSE_FLOAT
else
  PREC_TYPE = .double
endif

ifeq ($(FSANITIZE),TRUE)
  FSAN_TYPE = .san
  CXXFLAGS += -fsanitize=address,undefined
  LDFLAGS += -fsanitize=address,undefined
else
  FSAN_TYPE =
endif

BUILD_SUFFIX = $(BUILD_TYPE)$(CXX_TYPE)$(PREC_TYPE)$(FSAN_TYPE)
BUILD_DIR = tmp_build_dir/$(BUILD_SUFFIX)
EXE = shock$(BUILD_SUFFIX).ex

SRCS = main.cpp Hydro.cpp Riemann.cpp
OBJS = $(addprefix $(BUILD_DIR)/,$(SRCS:.cpp=.o))
DEPS = $(OBJS:.o=.d)

LDLIBS += -lm

default: $(EXE)

-include $(DEPS)

$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(DEPFLAGS) -c -o $@ $<

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

.PHONY: clean realclean

clean:
	$(RM) -r tmp_build_dir

realclean: clean
	$(RM) shock.*.ex
