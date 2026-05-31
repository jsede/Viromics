/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  darkMode: 'class', // support class-based dark mode
  theme: {
    extend: {
      fontFamily: {
        sans: ['Inter', 'sans-serif'],
        outfit: ['Inter', 'sans-serif'], // Map outfit to Inter to easily replace the font
        mono: ['Fira Code', 'Courier New', 'monospace'],
      },
      colors: {
        dark: {
          bg: '#F8FAFC',     // Light slate background
          card: '#FFFFFF',   // White panel cards
          border: '#E2E8F0', // Silver border
          hover: '#F1F5F9',  // Light hover
        },
        brand: {
          teal: '#DC2626',    // Lab Crimson Red
          violet: '#7C3AED',  // Secondary Violet
          emerald: '#10B981', // Clean Emerald
          rose: '#EF4444',    // Diagnostic Red
          indigo: '#4F46E5',  // Structural Indigo
          amber: '#F59E0B',   // Warn Amber
        }
      },
      boxShadow: {
        'glow-teal': '0 0 15px rgba(220, 38, 38, 0.15)', // Crimson glow
        'glow-violet': '0 0 15px rgba(124, 58, 237, 0.12)',
        'glow-rose': '0 0 15px rgba(239, 68, 68, 0.12)',
      }
    },
  },
  plugins: [],
}
