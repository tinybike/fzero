"use strict";

var test = require("tape");
var Decimal = require("decimal.js");
var fzero = require("../");

function cos(x) {
    return new Decimal(Math.cos(x.toNumber()));
}
function exp(x) {
    return new Decimal(Math.exp(-x.toNumber()) - 2);
}

test("fzero", function (t) {
    t.plan(2);
    t.equal(fzero(cos, 0, 3).x.toFixed(8), (Math.PI / 2).toFixed(8), "fzero(cos, 0, 3) == Pi/2");
    t.equal(fzero(exp, -10, 10).x.toFixed(8), "-0.69314718", "fzero(exp(-x) - 2, -10, 10) == -0.69314718");
    t.end();
});
