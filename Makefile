GP=$(HOME)/PARI/pari/gp

test:
	@$(GP) --test -q < test.in > /tmp/out
	@diff test.out /tmp/out
