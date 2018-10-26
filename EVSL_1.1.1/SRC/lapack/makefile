include ../../makefile.in

MUTE=@

# Object files
LAPACK_SRC = \
disnan.c    dlamch.c  dlar1v.c  dlarnv.c  dlarre.c  dlarrv.c  dlasq2.c  dlassq.c  dorg2r.c  dsteqr.c  dsytrd.c  zgetf2.c \
dlae2.c     dlaneg.c  dlarfb.c  dlarra.c  dlarrf.c  dlartg.c  dlasq5.c  dlatrd.c  dorgql.c  dsterf.c  ieeeck.c  zgetrf.c \
dlaebz.c    dlanst.c  dlarf.c   dlarrb.c  dlarrj.c  dlaruv.c  dlasq6.c  dlazq3.c  dorgqr.c  dstev.c   ilaenv.c  zgetrs.c \
dlaev2.c    dlansy.c  dlarfg.c  dlarrc.c  dlarrk.c  dlascl.c  dlasr.c   dlazq4.c  dorgtr.c  dsyev.c   iparmq.c  zlaswp.c \
dlaisnan.c  dlapy2.c  dlarft.c  dlarrd.c  dlarrr.c  dlaset.c  dlasrt.c  dorg2l.c  dstemr.c  dsytd2.c  zgesv.c

OBJS = ${LAPACK_SRC:.c=.o}

all: $(OBJS)

INCLUDES = -I. -I../../INC

clean:
	rm -rf *.o

distclean: clean
