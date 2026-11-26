/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      fontFamily: {
        sans: ["'Helvetica Neue'", 'Arial', '-apple-system', 'BlinkMacSystemFont', 'sans-serif'],
        mono: ["'SFMono-Regular'", 'Menlo', 'Consolas', "'Liberation Mono'", 'ui-monospace', 'monospace'],
      },
      colors: {
        surface: {
          canvas: '#F8F9FA',
          plate: '#FFFFFF',
          subtle: '#F1F3F5',
        },
        hairline: {
          DEFAULT: '#E2E6EA',
          strong: '#CBD2D9',
        },
        ink: {
          primary: '#111827',
          secondary: '#4B5563',
          muted: '#6B7280',
        },
        petrol: {
          DEFAULT: '#0F5C6E',
          hover: '#0B4654',
          tint: '#EAF3F5',
          border: '#B9D6DD',
        },
        qc: {
          pass: '#16654A',
          'pass-bg': '#ECF7F2',
          'pass-border': '#B5DFCE',
          'pass-bar': '#2A8C6A',
          warn: '#8C5308',
          'warn-bg': '#FDF6E7',
          'warn-border': '#EBD29B',
          'warn-bar': '#CC8400',
          fail: '#A62C22',
          'fail-bg': '#FCEEEB',
          'fail-border': '#EBBDB8',
          'fail-bar': '#C0392B',
        },
        bio: {
          missense: '#B43B2E',
          indel: '#9A5B06',
          synonymous: '#5C6773',
        },
      },
      boxShadow: {
        tooltip: '0 4px 12px rgba(17, 24, 39, 0.08)',
      },
    },
  },
  plugins: [],
}
