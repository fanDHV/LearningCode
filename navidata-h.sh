#! /usr/local/bin/tcsh -f
Navigator -h 25.off 25.in >& data/25.h.err
gprof Navigator gmon.out >& data/25.h.gprof.out

Navigator -h 50.off 50.in >& data/50.h.err
gprof Navigator gmon.out >& data/50.h.gprof.out

Navigator -h 100.off 100.in >& data/100.h.err
gprof Navigator gmon.out >& data/100.h.gprof.out

Navigator -h cone.off cone.in >& data/cone.h.err
gprof Navigator gmon.out >& data/cone.h.gprof.out

Navigator -h cone2.off cone2.in >& data/cone2.h.err
gprof Navigator gmon.out >& data/cone2.h.gprof.out

Navigator -h cone3.off cone3.in >& data/cone3.h.err
gprof Navigator gmon.out >& data/cone3.h.gprof.out

Navigator -h cone4.off cone4.in >& data/cone4.h.err
gprof Navigator gmon.out >& data/cone4.h.gprof.out

Navigator -h face19996.off face19996.in >& data/face19996.h.err
gprof Navigator gmon.out >& data/face19996.h.gprof.out

Navigator -h face49996.off face49996.in >& data/face49996.h.err
gprof Navigator gmon.out >& data/face49996.h.gprof.out

Navigator -h face5996.off face5996.in >& data/face5996.h.err
gprof Navigator gmon.out >& data/face5996.h.gprof.out

Navigator -h face7996.off face7996.in >& data/face7996.h.err
gprof Navigator gmon.out >& data/face7996.h.gprof.out

Navigator -h face99996.off face99996.in >& data/face99996.h.err
gprof Navigator gmon.out >& data/face99996.h.gprof.out

Navigator -h spiral1.off spiral1.in >& data/spiral1.h.err
gprof Navigator gmon.out >& data/spiral1.h.gprof.out
