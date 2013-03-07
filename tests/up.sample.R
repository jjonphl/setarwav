d4.s <- c(0.34150635, 0.59150635, 0.15849365, -0.09150635)

stopifnot(all.equal(c(d4.s[1],0,d4.s[2],0,d4.s[3],0,d4.s[4]),
                    up.sample(d4.s, 1)))
stopifnot(all.equal(c(d4.s[1],0,0,d4.s[2],0,0,d4.s[3],0,0,d4.s[4]),
                    up.sample(d4.s, 2)))
stopifnot(all.equal(c(d4.s[1],0,0,0,d4.s[2],0,0,0,d4.s[3],0,0,0,d4.s[4]),
                    up.sample(d4.s, 3)))
stopifnot(all.equal(c(d4.s[1],0,0,0,0,d4.s[2],0,0,0,0,d4.s[3],0,0,0,0,d4.s[4]),
                    up.sample(d4.s, 4)))
