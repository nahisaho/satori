import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    globals: true,
    testTimeout: 30000,
    include: ['tests/unit/**/*.test.js', 'tests/integration/**/*.test.js', 'tests/validation/**/*.test.js'],
    coverage: {
      provider: 'v8',
      include: ['bin/**/*.js'],
      reporter: ['text', 'json-summary'],
    },
  },
});
