import sys
import getopt

def run(x,y=5,z=True):
    print 'x = ',x
    print 'y = ',y
    print 'z = ',z

def run2(x,*argv):
    opts, args = getopt.getopt(argv,'abc:d:')
    print 'x = ',x
    print 'ops = ',opts
    print 'args = ',args
    print 'argv = ',argv

def run3(*argv):
    opts, args = getopt.getopt(argv,'abc:d:')
    print 'ops = ',opts
    print 'args = ',args
    print 'argv = ',argv


if __name__=='__main__':
    print sys.argv[1:]
    print tuple(sys.argv[1:])
    print 'run  *******'
    run(sys.argv[1:])
    print 'run  *******'
    run(*tuple(sys.argv[1:]))
    print 'run2 *******'
    run2(sys.argv[1:])
    print 'run3 *******'
    run3(*tuple(sys.argv[1:]))
    print'*'*30
