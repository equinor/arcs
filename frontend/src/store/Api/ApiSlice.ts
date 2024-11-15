import { createApi } from "@reduxjs/toolkit/query/react";
import { baseQueryWithReauth } from "./BaseQueryWithReauth";

export const apiSlice = createApi({
    reducerPath: "api",
    tagTypes: [
        "AvailableProjects",
        "Graph",
        "Components",
        "InputObject",
        "OperationalScenarioResults",
        "ProjectScenarios",
    ],
    baseQuery: baseQueryWithReauth,
    endpoints: () => ({}),
});
