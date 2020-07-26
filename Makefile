CC        := icc
CFLAGS    := -O3
DEPEND    := -MMD
LIBS      := -L$(HOME)/.local/lib -lfftw3
INCLUDES  := -I$(HOME)/.local/include -Iinclude
SRCSDIR   := src
OBJSDIR   := obj
SRCS      := $(foreach dir, $(shell find $(SRCSDIR) -type d), $(wildcard $(dir)/*.c))
OBJS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.o)))
DEPS      := $(addprefix $(OBJSDIR)/, $(subst $(SRCSDIR)/,,$(SRCS:.c=.d)))
TARGET    := a.out


$(TARGET): $(OBJS)
		$(CC) $(CFLAGS) $(DEPEND) -o $@ $^ $(LIBS) -lm

$(OBJSDIR)/%.o: $(SRCSDIR)/%.c
		@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
		$(CC) $(CFLAGS) $(DEPEND) $(LIBS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

.PHONY : clean
clean:
		$(RM) -r $(OBJSDIR)/* $(TARGET)

-include $(DEPS)

