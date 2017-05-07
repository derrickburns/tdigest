package(default_visibility = ["//visibility:public"])

cc_test(
    name = "tdigest_test",
    srcs = [
        "TDigestTest.cpp",
    ],
    size = "small",
    deps = [
        ":tdigest",
        "//external:gtest",
        "//external:glog",
    ],
    copts = [
        "-std=c++11",
    ],
)

cc_library(
    name = "tdigest",
    srcs = [
    ],
    hdrs = [
        "TDigest.h",
    ],
    deps = [
        "//external:glog",
        "//external:gtest",
        "//infra/serializer:serializer"
    ],
    copts = [
        "-std=c++11",
    ],
)


