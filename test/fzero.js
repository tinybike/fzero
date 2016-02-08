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

test("fzero", function (t) {
    t.plan(10);
    t.equal(fzero(cos, 0, 3).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3) == Pi/2");
    t.equal(fzero(cos, 0, 3, {eps: "0.05"}).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3, {eps: '0.05'}) == Pi/2");
    t.equal(fzero(cos, 0, 3, {eps: "0.05", mu: "0.25"}).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3, {eps: '0.05', mu: '0.25'}) == Pi/2");
    t.equal(fzero(cos, 0, 3, {maxiter: 50}).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3, {maxiter: 25}) == Pi/2");
    t.equal(fzero(cos, 0, 3, {maxiter: 50, maxfev: 45}).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3, {maxiter: 25, maxfev: 20}) == Pi/2");
    t.equal(fzero(exp, -10, 10).solution.toFixed(8), "-0.69314718", "fzero(exp(-x) - 2, -10, 10) == -0.69314718");
    t.equal(fzero(log, 0.01, 2).solution.toFixed(8), "1.00000000", "fzero(log, 0.01, 2) == 1.00000000");
    t.equal(fzero(hill, -1, 1).solution.toFixed(8), "-0.00000000", "fzero(-(x - 1)^2 + 1, -1, 1) == 0.00000000");
    t.equal(fzero(hill, 1, 3).solution.toFixed(8), "2.00000000", "fzero(-(x - 1)^2 + 1, 1, 3) == 2.00000000");
    t.equal(fzero(logistic, -10, 10).solution.toFixed(8), "-0.00000000", "fzero(1/(1 + e^(-x)) - 1/2, -10, 10) == 0.00000000");
    t.end();
});
