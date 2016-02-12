/**
 * Find a zero of a univariate function.
 *
 * This is essentially the ACM algorithm 748: Enclosing Zeros of
 * Continuous Functions due to Alefeld, Potra and Shi, ACM Transactions
 * on Mathematical Software, Vol. 21, No. 3, September 1995. Although
 * the workflow should be the same, the structure of the algorithm has
 * been transformed non-trivially; instead of the authors' approach of
 * sequentially calling building blocks subprograms we implement here a
 * FSM version using one interior point determination and one bracketing
 * per iteration, thus reducing the number of temporary variables and
 * simplifying the algorithm structure. Further, this approach reduces
 * the need for external functions and error handling. The algorithm has
 * also been slightly modified.
 *
 * (This is a port of Octave's fzero implementation, by Jaroslav Hajek.)
 *
 * @author Jack Peterson (jack@tinybike.net)
 */

"use strict";

var Decimal = require("decimal.js");

function unique(arr) {
    var u = {}, a = [];
    for (var i = 0, l = arr.length; i < l; ++i) {
        if (u.hasOwnProperty(arr[i])) continue;
        a.push(arr[i]);
        u[arr[i]] = 1;
    }
    return a;
}

function toDecimal(x) {
    if (x && x.constructor !== Decimal) {
        if (x.toFixed && x.toFixed.constructor === Function) x = x.toFixed();
        x = new Decimal(x);
    }
    return x;
}

module.exports = function (f, bounds, options) {
    options = options || {};
    var mu = toDecimal(options.mu) || new Decimal("0.5");
    var eps = toDecimal(options.eps) || new Decimal("0.001");
    var tolx = toDecimal(options.tolx) || new Decimal(0);
    var maxiter = options.maxiter || 100;
    var maxfev = options.maxfev || maxiter;
    var verbose = options.verbose;

    // The default exit flag if exceeded number of iterations.
    var code = 0;
    var niter = 0;
    var nfev = 0;

    var x = new Decimal(NaN);
    var fval = new Decimal(NaN);
    var a = new Decimal(NaN);
    var fa = new Decimal(NaN);
    var b = new Decimal(NaN);
    var fb = new Decimal(NaN);

    // Prepare...
    if (bounds === null || bounds === undefined) {
        throw new Error("Initial guess required");
    }
    if (bounds.constructor === Array && bounds.length) {
        a = new Decimal(bounds[0].toString());
    } else {
        a = new Decimal(bounds.toString());
    }
    fa = toDecimal(f(a.toString()));
    nfev = 1;
    if (bounds.constructor === Array && bounds.length > 1) {
        b = new Decimal(bounds[1].toString());
        fb = toDecimal(f(b.toString()));
        nfev += 1;
    } else {
        // Try to get b.
        if (verbose) {
            console.log("Search for an interval around", a.toString(), "containing a sign change:");
            console.log("count\ta\t\tf(a)\t\t\t\tb\t\tf(b)");
        }
        var bracketed, blist;
        var aa = (a.eq(new Decimal(0))) ? new Decimal(1) : a;
        var tries = 0;
        do {
            blist = [
                aa.times(new Decimal("0.9")),
                aa.times(new Decimal("1.1")),
                aa.minus(new Decimal(1)),
                aa.plus(new Decimal(1)),
                aa.times(new Decimal("0.5")),
                aa.times(new Decimal("1.5")),
                aa.neg(),
                aa.times(new Decimal(2)),
                aa.times(new Decimal(10)).neg(),
                aa.times(new Decimal(10))
            ];
            for (var j = 0, len = blist.length; j < len; ++j) {
                b = blist[j];
                fb = toDecimal(f(b.toString()));
                if (verbose) {
                    console.log(nfev + "\t\t" + aa + "\t\t" + fa + "\t\t" + b + "\t\t" + fb);
                }
                nfev += 1;
                if (fa.s * fb.s <= 0) {
                    a = aa;
                    bracketed = true;
                    break;
                }
            }
            if (!options.randomize) break;
            aa = aa.times(new Decimal(Math.random().toString())).plus(new Decimal((Math.random() - 0.5).toString()));
            fa = toDecimal(f(aa.toString()));
        } while (!bracketed && ++tries < maxiter);
    }

    var u, fu;
    if (b.lt(a)) {
        u = a;
        a = b;
        b = u;
        fu = fa;
        fa = fb;
        fb = fu;
    }
    if (fa.s * fb.s > 0) {
        throw new Error("Invalid initial bracketing");
    }

    var slope0 = fb.minus(fa).dividedBy(b.minus(a));
    if (fa.eq(new Decimal(0))) {
        b = a;
        fb = fa;
    } else if (fb.eq(new Decimal(0))) {
        a = b;
        fa = fb;
    }
    var itype = 1;
    if (fa.abs().lt(fb.abs())) {
        u = a;
        fu = fa;
    } else {
        u = b;
        fu = fb;
    }
    var d = u;
    var e = u;
    var fd = fu;
    var fe = fu;
    var mba = mu.times(b.minus(a));
    var c, fc, df, procedure;
    if (verbose) {
        console.log("Search for a zero in the interval [" + a.toString() + ", " + b.toString() + "]:");
        console.log("count\tx\t\t\tf(x)\t\t\tprocedure");
        console.log(nfev + "\t\t" + a.toFixed(18) + "\t" + fa.toFixed(18) + "\tinitial lower");
        console.log(nfev + "\t\t" + b.toFixed(18) + "\t" + fb.toFixed(18) + "\tinitial upper");
    }
    while (niter < maxiter && nfev < maxfev) {
        switch (itype) {
        case 1:
            // The initial test.
            if (b.minus(a).lte(u.abs().times(eps).times(new Decimal(2)).plus(tolx).times(new Decimal(2)))) {
                x = u;
                fval = fu;
                code = 1;
            } else {
                if (fa.abs().lte(fb.abs().times(new Decimal(1000))) && (fb.abs().lte(fa.abs().times(new Decimal(1000))))) {
                    // Secant step.
                    c = u.minus(a.minus(b).dividedBy(fa.minus(fb)).times(fu));
                    procedure = "secant";
                } else {
                    // Bisection step.
                    c = a.plus(b).dividedBy(new Decimal(2));
                    procedure = "bisection";
                }
                d = u;
                df = fu;
                itype = 5;
            }
            break;
        case 2:
        case 3:
            var l = unique([fa, fb, fd, fe]).length;
            if (l === 4) {
                // Inverse cubic interpolation.
                var q11 = d.minus(e).times(fd).dividedBy(fe.minus(fd));
                var q21 = b.minus(d).times(fb).dividedBy(fd.minus(fb));
                var q31 = a.minus(b).times(fa).dividedBy(fb.minus(fa));
                var d21 = b.minus(d).times(fd).dividedBy(fd.minus(fb));
                var d31 = a.minus(b).times(fb).dividedBy(fb.minus(fa));
                var q22 = d21.minus(q11).times(fb).dividedBy(fe.minus(fb));
                var q32 = d31.minus(q21).times(fa).dividedBy(fd.minus(fa));
                var d32 = d31.minus(q21).times(fd).dividedBy(fd.minus(fa));
                var q33 = d32.minus(q22).times(fa).dividedBy(fe.minus(fa));
                c = a.plus(q31).plus(q32).plus(q33);
                procedure = "interpolation (cubic)";
            }
            if (l < 4 || c.minus(a).s * c.minus(b).s < 0) {
                // Quadratic interpolation + Newton.
                var a0 = fa;
                var a1 = fb.minus(fa).dividedBy(b.minus(a));
                var a2 = (fd.minus(fb).dividedBy(d.minus(b)).minus(a1)).dividedBy(d.minus(a));
                // Modification 1
                c = a.minus(a0.dividedBy(a1));
                if (!a2.eq(new Decimal(0))) {
                    c = a.minus(a0.dividedBy(a1));
                    var pc, pdc;
                    for (var k = 0; k < itype; ++k) {
                        pc = a0.plus(a1.plus(a2.times(c.minus(b)).times(c.minus(a))));
                        pdc = a1.plus(a2.times(c.times(new Decimal(2)).minus(a).minus(b)));
                        if (pdc.eq(new Decimal(0))) {
                            c = a.minus(a0.dividedBy(a1));
                        } else {
                            c = c.minus(pc.dividedBy(pdc));
                        }
                    }
                }
                procedure = "interpolation (quadratic)";
            }
            itype++;
            break;
        case 4:
            // Double secant step.
            c = u.minus(b.minus(a).times(new Decimal(2)).dividedBy(fb.minus(fa)).times(fu));
            if (c.minus(u).abs().gt(b.minus(a).dividedBy(new Decimal(2)))) {
                c = b.plus(a).dividedBy(new Decimal(2));
            }
            itype = 5;
            break;
        case 5:
            // Bisection step.
            c = b.plus(a).dividedBy(new Decimal(2));
            procedure = "bisection";
            itype = 2;
        }

        // Don't let c come too close to a or b.
        var delta = u.abs().times(eps).plus(tolx).times(new Decimal("0.14"));
        if (b.minus(a).lte(delta.times(new Decimal(2)))) {
            c = a.plus(b).dividedBy(new Decimal(2));
        } else {
            var s;
            if (b.minus(delta).lt(c)) {
                s = b.minus(delta);
            } else {
                s = c;
            }
            if (a.plus(delta).gt(s)) {
                c = a.plus(delta);
            } else {
                c = s;
            }
        }

        // Calculate new point.
        x = c;
        fval = toDecimal(f(c.toString()));
        fc = fval;
        niter++;
        nfev++;
        if (verbose) {
            console.log(nfev + "\t\t" + x.toFixed(18) + "\t" + fval.toFixed(18) + "\t" + procedure);
        }

        // Modification 2: skip inverse cubic interpolation if non-monotonicity
        // is detected.
        if (fc.minus(fa).s * fc.minus(fb).s >= 0) {
            // The new point broke monotonicity; disable inverse cubic.
            fe = fc;
        } else {
            e = d;
            fe = fd;
        }

        // Bracketing
        if (fa.s * fc.s < 0) {
            d = b;
            fd = fb;
            b = c;
            fb = fc;
        } else if (fb.s * fc.s < 0) {
            d = a;
            fd = fa;
            a = c;
            fa = fc;
        } else if (fc.eq(new Decimal(0))) {
            b = c;
            a = b;
            fb = fc;
            fa = fb;
            code = 1;
        } else {
            // This should never happen.
            throw new Error("Zero point is not bracketed");
        }

        if (fa.abs().lt(fb.abs())) {
            u = a;
            fu = fa;
        } else {
            u = b;
            fu = fb;
        }
        if (b.minus(a).lte(u.abs().times(eps).times(new Decimal(2)).plus(tolx))) {
            code = 1;
        }

        // Skip bisection step if successful reduction.
        if (itype === 5 && b.minus(a).lte(mba)) {
            itype = 2;
        }
        if (itype === 2) {
            mba = mu.times(b.minus(a));
        }

    } // while

    // Check solution for a singularity by examining slope.
    if (code === 1) {
        var m;
        if (new Decimal(1e6).gt(new Decimal("0.5").dividedBy(eps.plus(tolx)))) {
            m = new Decimal(1e6);
        } else {
            m = new Decimal("0.5").dividedBy(eps.plus(tolx));
        }
        if (!b.minus(a).eq(new Decimal(0)) && fb.minus(fa).dividedBy(b.minus(a)).dividedBy(slope0).gt(m)) {
            code = -5;
        }
    }

    return {
        solution: x.toFixed(),
        fval: fval.toString(),
        code: code,
        diagnostic: {
            iterations: niter,
            functionEvals: nfev,
            bracketx: [a.toFixed(), b.toFixed()],
            brackety: [fa.toString(), fb.toString()]
        }
    };
};
