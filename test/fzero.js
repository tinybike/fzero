"use strict";

var test = require("tape");
var Decimal = require("decimal.js");
var fzero = require("../");

var ONE = new Decimal(1);
var TWO = new Decimal(2);
var HALF = ONE.dividedBy(TWO);

function cos(x) {
    return Math.cos(Number(x)).toString();
}
function exp(x) {
    return (Math.exp(-Number(x)) - 2).toString();
}
function log(x) {
    return Math.log(Number(x)).toString();
}
function hill(x) {
    return new Decimal(x).minus(ONE).toPower(TWO).neg().plus(ONE);
}
function logistic(x) {
    return ONE.dividedBy(ONE.plus(new Decimal(x).neg().exp())).minus(HALF);
}
function singular(x) {
    return ONE.dividedBy(new Decimal(x).minus(ONE));
}
function lslmsr(n, q, i, a, xi) {
    n = new Decimal(n);
    var numOutcomes = q.length;
    var qj = new Array(numOutcomes);
    var sum_q = new Decimal(0);
    for (var j = 0; j < numOutcomes; ++j) {
        qj[j] = q[j];
        sum_q = sum_q.plus(q[j]);
    }
    qj.splice(i, 1);
    var q_plus_n = n.plus(sum_q);
    var b = a.times(q_plus_n);
    var exp_qi = q[i].plus(n).dividedBy(b).exp();
    var exp_qj = new Array(numOutcomes);
    var sum_qj = new Decimal(0);
    var sum_exp_qj = new Decimal(0);
    var sum_qj_x_expqj = new Decimal(0);
    for (j = 0; j < numOutcomes - 1; ++j) {
        sum_qj = sum_qj.plus(qj[j]);
        exp_qj[j] = qj[j].dividedBy(b).exp();
        sum_exp_qj = sum_exp_qj.plus(exp_qj[j]);
        sum_qj_x_expqj = sum_qj_x_expqj.plus(q[j].times(exp_qj[j]));
    }
    return a.times(q[i].plus(n).dividedBy(b).exp().plus(sum_exp_qj).ln()).plus(
        exp_qi.times(sum_qj).minus(sum_qj_x_expqj).dividedBy(
            q_plus_n.times(exp_qi.plus(sum_exp_qj))
        ).minus(xi)
    );
}


test("fzero", function (t) {
    var trial = function (tr) {
        var code = (tr.code === undefined) ? 1 : tr.code;
        var options = (tr.options === undefined) ? "" : ", " + JSON.stringify(tr.options);
        var actual = fzero(tr.f, tr.bounds, tr.options);
        t.equal(Number(actual.solution).toFixed(8), tr.expected, "fzero(" + tr.label + ", " + JSON.stringify(tr.bounds) + options + ") == " + tr.expected);
        t.equal(actual.code, code, "exit code == " + code);
    };
    trial({
        f: cos,
        label: "cos",
        bounds: [0, 3],
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: 1.0,
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [1.0],
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [1.0],
        options: {randomize: true},
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [0, 3],
        options: {eps: "0.05"},
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [0, 3],
        options: {eps: "0.05", mu: "0.25"},
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [0, 3],
        options: {maxiter: 50},
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: cos,
        label: "cos",
        bounds: [0, 3],
        options: {maxiter: 50, maxfev: 45},
        expected: (Math.PI / 2).toFixed(8)
    });
    trial({
        f: exp,
        label: "exp",
        bounds: [-10, 10],
        expected: "-0.69314718"
    });
    trial({
        f: log,
        label: "log",
        bounds: [0.01, 2],
        expected: "1.00000000"
    });
    trial({
        f: hill,
        label: "-(x - 1)^2 + 1",
        bounds: [-1, 1],
        expected: "-0.00000000",
        code: 0
    });
    trial({
        f: hill,
        label: "-(x - 1)^2 + 1",
        bounds: [1, 3],
        expected: "2.00000000"
    });
    trial({
        f: hill,
        label: "-(x - 1)^2 + 1",
        bounds: [2, 3],
        expected: "2.00000000"
    });
    trial({
        f: hill,
        label: "-(x - 1)^2 + 1",
        bounds: [1, 2.00000000000001],
        expected: "2.00000000"
    });
    trial({
        f: logistic,
        label: "1/(1 + e^(-x)) - 1/2",
        bounds: [-10, 10],
        expected: "-0.00000000",
        code: -5
    });
    trial({
        f: singular,
        label: "1/(x - 1)",
        bounds: [-10, 10],
        expected: "1.00000000",
        code: -5
    });
    trial({
        f: function (n) {
            var q = [new Decimal(10),      // outcome 1 shares
                     new Decimal(10),      // outcome 2 shares
                     new Decimal(10),      // outcome 3 shares
                     new Decimal(10),      // outcome 4 shares
                     new Decimal(10)];     // outcome 5 shares
            var i = 0;                     // array index of outcome to trade
            var a = new Decimal("0.0079"); // LS-LMSR alpha
            var xi = new Decimal("0.3");   // price cap
            return lslmsr(n, q, i, a, xi);
        },
        label: "lslmsr(n, [10, 10, 10, 10, 10], 0, '0.0079', '0.3')",
        bounds: [1e-12, 1000],
        expected: "0.18973664"
    });
    trial({
        f: function (n) {
            var q = [new Decimal("659.90262467263840222037"),
                     new Decimal("666.57262467263840222039")];
            var i = 0;
            var a = new Decimal("0.00790000000000000001");
            var xi = new Decimal("0.5");
            return lslmsr(n, q, i, a, xi);
        },
        label: "lslmsr(n, ['659.90262467263840222037', '666.57262467263840222039'], 0, '0.00790000000000000001', '0.5')",
        bounds: [1e-12, 1000],
        expected: "6.33231266"
    })
    trial({
        f: function (n) {
            var q = [new Decimal(10),      // outcome 1 shares
                     new Decimal(10),      // outcome 2 shares
                     new Decimal(10),      // outcome 3 shares
                     new Decimal(10),      // outcome 4 shares
                     new Decimal(10)];     // outcome 5 shares
            var i = 0;                     // array index of outcome to trade
            var a = new Decimal("0.0079"); // LS-LMSR alpha
            var xi = new Decimal("0.3");   // price cap
            var ans = lslmsr(n, q, i, a, xi);
            return ans;
        },
        label: "lslmsr(n, [10, 10, 10, 10, 10], 0, '0.0079', '0.3')",
        bounds: [1e-12, 1000],
        expected: "0.18973664"
    });
    trial({
        f: function (n) {
            var q = [new Decimal("71.43667960674091200687"),
                     new Decimal("72.6286796067409120068")];
            var xi = new Decimal("0.6");
            var a = new Decimal("0.00790000000000000001");
            var i = 0;
            var ans = lslmsr(n, q, i, a, xi);
            return ans;
        },
        label: "lslmsr(n, ['71.43667960674091200687', '72.6286796067409120068'], 0, '0.00790000000000000001', '0.6')",
        bounds: [1e-12, 1000],
        options: {verbose: true},
        expected: "1.61713310"
    });
    trial({
        f: function (n) {
            var q = [new Decimal("71.43667960674091200687"),
                     new Decimal("72.6286796067409120068")];
            var xi = new Decimal("0.6");
            var a = new Decimal("0.00790000000000000001");
            var i = 0;
            var ans = lslmsr(n, q, i, a, xi);
            return ans;
        },
        label: "lslmsr(n, ['71.43667960674091200687', '72.6286796067409120068'], 0, '0.00790000000000000001', '0.6')",
        bounds: [-0.78, 1.78],
        options: {tolx: 1e-12, maxiter: 1000},
        expected: "1.61713310"
    });
    t.throws(function () { fzero(log, [2, 3]); }, /Invalid initial bracketing/, "fzero(log, 2, 3) throws Error('Invalid initial bracketing')");
    t.throws(function () { fzero(log, [-1, 0]); }, /Zero point is not bracketed/, "fzero(log, -1, 0) throws Error('Zero point not bracketed')");
    t.throws(function () { fzero(exp); }, /Initial guess required/, "fzero(exp)");
    t.end();
});
