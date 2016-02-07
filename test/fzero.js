"use strict";

var test = require("tape");
var Decimal = require("decimal.js");
var fzero = require("../");

function cos(x) {
    return Math.cos(Number(x)).toString();
}
function exp(x) {
    return (Math.exp(-Number(x)) - 2).toString();
}

test("fzero", function (t) {
    t.plan(2);
    t.equal(fzero(cos, 0, 3).solution.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3) == Pi/2");
    t.equal(fzero(exp, -10, 10).solution.toFixed(8), "-0.69314718", "fzero(exp(-x) - 2, -10, 10) == -0.69314718");
    t.end();
});
