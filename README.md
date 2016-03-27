fzero
=====

[![Build Status](https://travis-ci.org/tinybike/fzero.svg)](https://travis-ci.org/tinybike/fzero)
[![Coverage Status](https://coveralls.io/repos/tinybike/fzero/badge.svg?branch=master&service=github)](https://coveralls.io/github/tinybike/fzero?branch=master)
[![npm version](https://badge.fury.io/js/fzero.svg)](https://badge.fury.io/js/fzero)

![Falcon Punch](/falcon.jpg?raw=true "Falcon Punch")

Find a zero of a univariate function.  This is a JavaScript port of [Jaroslav Hajek](highegg@gmail.com)'s [fzero implementation](https://fossies.org/dox/octave-4.0.0/fzero_8m_source.html) in [Octave](https://www.gnu.org/software/octave/).

This is essentially the ACM algorithm 748: Enclosing Zeros of Continuous Functions due to Alefeld, Potra and Shi, ACM Transactions on Mathematical Software, Vol. 21, No. 3, September 1995.  Although the workflow should be the same, the structure of the algorithm has been transformed non-trivially; instead of the authors' approach of sequentially calling building blocks subprograms we implement here a FSM version using one interior point determination and one bracketing per iteration, thus reducing the number of temporary variables and simplifying the algorithm structure.  Further, this approach reduces the need for external functions and error handling.  The algorithm has also been slightly modified.

Usage
-----
```
$ npm install fzero
```
To use fzero in Node.js, simply require it:
```javascript
var fzero = require("fzero");
```
A minified, browserified file `dist/fzero.min.js` is included for use in the browser.  Including this file attaches a `fzero` object to `window`:
```html
<script src="dist/fzero.min.js" type="text/javascript"></script>
```
`fzero` is a function that takes 3 required arguments: the function to find the zero of, a lower-bound for the zero, and an upper-bound for the zero.  A fourth optional argument can be used to adjust fzero's settings; for example, `maxiter: 50` sets the maximum number of iterations to 50, and `verbose: true` prints details of each iteration.  (See tests for details.)
```javascript
var myFunction = function (x) { return Math.cos(Number(x)).toString(); }
var lowerBound = 0;
var upperBound = 3;
var options = {maxiter: 50};
var zero = fzero(myFunction, [lowerBound, upperBound], options);
```
If you don't know the exact bounds, you can simply pass `fzero` an initial guess instead:
```javascript
var initialGuess = 2;
var zero = fzero(myFunction, initialGuess, options);
```
`fzero` returns an object that has the following fields:

- `solution`: `x` such that `myFunction(x) = 0`
- `fval`: The numerical value of `myFunction` at `solution`.
- `code`: Exit flag which can have one of the following values.
    - `1`: The algorithm converged to a solution.
    - `0`: Maximum number of iterations or function evaluations has been reached.
    - `-1`: The algorithm has been terminated from user output function.
    - `-5`: The algorithm may have converged to a singular point.
- `diagnostic`: An object which has the following fields.
    - `iterations`: The number of iterations completed.
    - `functionEvals`: Number of times `myFunction` was evaluated.
    - `bracketx`: An array `[lower, upper]` with the ending lower and upper x-bounds.
    - `brackety`: An array `[f(lower), f(upper)]` with the ending lower and upper y-bounds.

Note: `fzero` uses [decimal.js](https://github.com/MikeMcl/decimal.js/) for arithmetic, so the input function should accept a string input (rather than a JS number).

Tests
-----
Unit tests are included in `test/`, and can be run using npm:
```
$ npm test
```
