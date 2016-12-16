CFLAGS = -Wall -Wextra -g
LDLIBS = -lm4ri

all: test_isd partial_col

test_isd: test_isd.c
	  $(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

partial_col: partial_col.c
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)
