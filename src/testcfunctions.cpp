#include <testthat.h>
#include "largeVis.h"
#include "alias.h"
#include "gradients.h"

// Initialize a unit test context. This is similar to how you
// might begin an R test file with 'context()', expect the
// associated context should be wrapped in braced.
context("alias tests") {

  test_that("alias int succeeds") {
  	AliasTable<int, double, double> testAlias = AliasTable<int, double, double>(100);
  	double weights[100];
  	for (int i = 0; i != 100; i++) weights[i] = i;
  	testAlias.initialize(weights);
  	testAlias.initRandom(1024);
  	expect_true(testAlias() == 71);
  	expect_true(testAlias() == 74);
  	expect_true(testAlias() == 70);
  	expect_true(testAlias() == 90);
  }

  test_that("alias long succeeds") {
  	AliasTable<long, double, double> testAlias = AliasTable<long, double, double>(100);
  	double weights[100];
  	for (int i = 0; i != 100; i++) weights[i] = i;
    testAlias.initialize(weights);
    testAlias.initRandom(1024);
    expect_true(testAlias() == 71);
    expect_true(testAlias() == 74);
    expect_true(testAlias() == 70);
    expect_true(testAlias() == 90);
  }

	test_that("alias long long succeeds") {
		AliasTable<long long, double, double> testAlias = AliasTable<long long, double, double>(100);
		double weights[100];
		for (int i = 0; i != 100; i++) weights[i] = i;
		testAlias.initialize(weights);
		testAlias.initRandom(1024);
		expect_true(testAlias() == 71);
		expect_true(testAlias() == 74);
		expect_true(testAlias() == 70);
		expect_true(testAlias() == 90);
}

	test_that("alias vertexidxtype succeeds") {
		AliasTable<vertexidxtype, double, double> testAlias = AliasTable<vertexidxtype, double, double>(100);
		double weights[100];
		for (int i = 0; i != 100; i++) weights[i] = i;
		testAlias.initialize(weights);
		testAlias.initRandom(1024);
		expect_true(testAlias() == 71);
		expect_true(testAlias() == 74);
		expect_true(testAlias() == 70);
		expect_true(testAlias() == 90);
	}

	test_that("alias edgeidxtype succeeds") {
		AliasTable<edgeidxtype, double, double> testAlias = AliasTable<edgeidxtype, double, double>(100);
		double weights[100];
		for (int i = 0; i != 100; i++) weights[i] = i;
		testAlias.initialize(weights);
		testAlias.initRandom(1024);
		expect_true(testAlias() == 71);
		expect_true(testAlias() == 74);
		expect_true(testAlias() == 70);
		expect_true(testAlias() == 90);
	}
}

context("gradient tests") {

	double x_i[2], y_i[2], holder[2];

	x_i[0] = -1; x_i[1] = 1;
	y_i[0] = -2; y_i[1] = 2;
	holder[0] = holder[1] = 0;

	AlphaGradient alphaGrad = AlphaGradient(2, 5, 2);
	test_that("alpha gradient positive") {
		alphaGrad.positiveGradient(x_i, y_i, holder);
		expect_true(holder[0] == -0.8);
		expect_true(holder[1] == - holder[0]);
	}
	test_that("alpha gradient negative") {
		alphaGrad.negativeGradient(x_i, y_i, holder);
		expect_true(holder[0] == 2);
		expect_true(holder[1] == - holder[0]);
	}

	AlphaOneGradient aOne = AlphaOneGradient(5, 2);
	test_that("alpha one gradient positive") {
		aOne.positiveGradient(x_i, y_i, holder);
		expect_true(pow(holder[0] - -0.66666667, 2) < .2);
		expect_true(holder[1] == - holder[0]);
	}
	test_that("alpha one gradient negative") {
		aOne.negativeGradient(x_i, y_i, holder);
		expect_true(pow(holder[0] - 1.5873015873, 2) < 0.1);
		expect_true(holder[1] == - holder[0]);
	}

	ExpGradient expG = ExpGradient(5, 2);
	test_that("exp gradient positive") {
		expG.positiveGradient(x_i, y_i, holder);
		expect_true(pow(holder[0] - -0.880797078, 2) < 0.1);
		expect_true(holder[1] == - holder[0]);
	}
	test_that("exp gradient negative") {
		expG.negativeGradient(x_i, y_i, holder);
		expect_true(pow(holder[0] -0.5960146101, 2) < 0.1);
		expect_true(holder[1] == - holder[0]);
	}
}
